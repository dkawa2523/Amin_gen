from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from pathlib import Path

from gas_screening_mvp.config_loader import load_config
from gas_screening_mvp.derivation.candidate_annotations import (
    amine_class_summary_rows,
    amine_summary_rows,
    candidate_annotation_row,
    candidate_breakdown_rows,
    candidate_classification_records,
)
from gas_screening_mvp.derivation.summary_builder import build_summary_row
from gas_screening_mvp.export.csv_exporter import export_csv
from gas_screening_mvp.export.excel_exporter import export_excel
from gas_screening_mvp.generation.registry import generate_candidates
from gas_screening_mvp.normalization.dedupe import dedupe_molecules
from gas_screening_mvp.normalization.rdkit_normalizer import normalize_all
from gas_screening_mvp.planning.fetch_planner import FetchDecision, FetchPlanner
from gas_screening_mvp.prefilter.structure_filter import CandidatePrefilter
from gas_screening_mvp.providers.amine_estimates import AmineLocalEstimateProvider
from gas_screening_mvp.providers.chemicals_local import ChemicalsLocalProvider
from gas_screening_mvp.providers.coolprop_local import CoolPropLocalProvider
from gas_screening_mvp.providers.curated_csv import CuratedCsvPropertyProvider
from gas_screening_mvp.providers.gwp_local import GwpLocalProvider
from gas_screening_mvp.providers.kinetics_probe import LocalKineticsProbeProvider
from gas_screening_mvp.providers.local_identity import LocalIdentityProvider
from gas_screening_mvp.providers.persistence import PersistenceScreener
from gas_screening_mvp.providers.pfas_local import PfasLocalClassifier
from gas_screening_mvp.providers.pubchem_pugrest import PubChemPugRestProvider
from gas_screening_mvp.providers.pubchem_pugview import PubChemPugViewProvider
from gas_screening_mvp.providers.reactivity import ReactivityClassifier
from gas_screening_mvp.providers.thermo_vapor_pressure import ThermoVaporPressureProvider
from gas_screening_mvp.selection.property_selector import PropertySelector
from gas_screening_mvp.storage.cache import SqliteApiCache


MODES = {"exploration", "enrichment", "refresh"}


class ScreeningPipeline:
    def __init__(self, config: dict):
        self.config = deepcopy(config)
        self.mode = _normalize_mode(self.config.get("mode"))
        self.config["mode"] = self.mode
        providers_cfg = self.config["providers"]
        self.remote_enabled = bool(providers_cfg.get("pubchem_enabled", False) or providers_cfg.get("pugview_enabled", False))
        self.pubchem_enabled = bool(providers_cfg.get("pubchem_enabled", False))
        self.pugview_enabled = bool(providers_cfg.get("pugview_enabled", False))
        self.max_api_candidates_per_run = int(providers_cfg.get("max_api_candidates_per_run", 25))

        cache_db = providers_cfg.get("cache_db") or "./gas_screening_cache.sqlite"
        self.cache = SqliteApiCache(cache_db)
        self.fetch_planner = FetchPlanner(self.cache, self.config.get("fetch_planning", {}))
        self.local_identity = LocalIdentityProvider()
        self.pubchem = PubChemPugRestProvider(
            cache=self.cache,
            enabled=self.pubchem_enabled and self.remote_enabled,
            rps=float(providers_cfg.get("pubchem_rps", 3.0)),
        )

        target_T = self.config["selection"].get("target_vapor_temperatures_K", [298.15, 313.15, 333.15])
        self.property_providers = [
            CuratedCsvPropertyProvider(providers_cfg.get("curated_properties_csv")),
            CoolPropLocalProvider(target_T),
            ChemicalsLocalProvider(),
            ThermoVaporPressureProvider(target_T),
            AmineLocalEstimateProvider(target_T),
            GwpLocalProvider(providers_cfg.get("gwp_csv")),
        ]
        self.pugview = PubChemPugViewProvider(
            cache=self.cache,
            enabled=self.pugview_enabled and self.remote_enabled,
            rps=float(providers_cfg.get("pugview_rps", 1.0)),
        )
        self.pfas = PfasLocalClassifier(providers_cfg.get("pfas_list_csv"))
        self.reactivity = ReactivityClassifier(providers_cfg.get("reactivity_csv"))
        self.kinetics = LocalKineticsProbeProvider(providers_cfg.get("kinetics_csv"))
        self.selector = PropertySelector(
            provider_priority=self.config["selection"].get("provider_priority"),
            target_units=self.config["selection"].get("target_units"),
        )
        self.stats = self._empty_stats()

    def run(self, input_csv: str | Path | None = None, dry_run: bool = False) -> dict[str, list[dict]]:
        candidates = generate_candidates(
            self.config["generation"],
            input_csv,
            cache=self.cache,
            mode=self.mode,
            remote_enabled=self.remote_enabled,
            dry_run=dry_run,
        )
        normalized = normalize_all(
            candidates,
            allowed_elements=set(self.config["prefilter"].get("allowed_elements", [])) or None,
        )
        unique = dedupe_molecules(normalized)

        prefilter = CandidatePrefilter(self.config["prefilter"])
        passed_results, rejected_results = prefilter.filter(unique)
        passed = [r.molecule for r in passed_results]

        self.stats.update({
            "generated_candidates": len(candidates),
            "normalized_candidates": len(normalized),
            "deduplicated_candidates": len(unique),
            "prefilter_passed": len(passed_results),
            "prefilter_rejected": len(rejected_results),
        })

        rejected_rows = [
            {
                "candidate_id": r.molecule.candidate_id,
                "input_name": r.molecule.input_name,
                "family": r.molecule.family,
                "score": r.api_priority_score,
                "reasons": "; ".join(r.reasons),
            }
            for r in rejected_results
        ]

        if dry_run:
            planned_rows = self._planned_api_rows(passed, apply_limit=True)
            self._record_planned_api_stats(planned_rows)
            self.stats["final_summary_rows"] = 0
            return self._sheets([], [], [], [], rejected_rows, planned_rows)

        summary_rows: list[dict] = []
        evidence_rows: list[dict] = []
        review_rows: list[dict] = []
        planned_rows: list[dict] = []
        remote_candidates_used = 0

        for mol in passed:
            chemical = self.local_identity.resolve(mol)[0]
            property_candidates = self._fetch_properties(chemical)
            selected = self._select_properties(chemical, property_candidates)
            decisions = self._plan_remote(mol, selected)

            if decisions and remote_candidates_used >= self.max_api_candidates_per_run:
                decisions = []
            if decisions:
                remote_candidates_used += 1
                planned_rows.extend(self._decision_rows(mol.candidate_id, decisions))

            remote_identity = self._resolve_remote_identity(mol, decisions)
            if remote_identity is not None:
                chemical = remote_identity
                property_candidates = self._fetch_properties(chemical)
                selected = self._select_properties(chemical, property_candidates)

            pfas_records = self.pfas.classify(chemical)
            persistence = PersistenceScreener({chemical.chemical_id: pfas_records}).classify(chemical)
            reactivity = (
                self.reactivity.classify(chemical)
                + candidate_classification_records(mol, chemical.chemical_id)
                + self._classify_pugview(chemical, decisions)
            )
            classifications = pfas_records + persistence + reactivity
            kinetics = self.kinetics.probe(chemical, self.config["reaction_probe_targets"])
            summary = build_summary_row(
                chemical,
                selected,
                classifications,
                kinetics,
                thresholds=self.config.get("supply_thresholds", {}),
            )
            summary.update(candidate_annotation_row(mol))
            summary_rows.append(summary)
            evidence_rows.extend(self._evidence_rows(chemical, property_candidates, selected, classifications, kinetics))
            if summary.get("review_required"):
                review_rows.append(self._review_row(summary))

        self.stats["remote_enrichment_candidates"] = remote_candidates_used
        self._record_planned_api_stats(planned_rows)
        self.stats["final_summary_rows"] = len(summary_rows)
        return self._sheets(summary_rows, evidence_rows, self._coverage_rows(summary_rows), review_rows, rejected_rows, planned_rows)

    def run_to_excel(self, output_xlsx: str | Path, input_csv: str | Path | None = None, dry_run: bool = False) -> dict[str, list[dict]]:
        sheets = self.run(input_csv=input_csv, dry_run=dry_run)
        export_excel(sheets, output_xlsx)
        return sheets

    def run_to_output(self, output_path: str | Path, input_csv: str | Path | None = None, dry_run: bool = False) -> dict[str, list[dict]]:
        sheets = self.run(input_csv=input_csv, dry_run=dry_run)
        output = Path(output_path)
        if output.suffix.lower() in {".xlsx", ".xlsm"}:
            export_excel(sheets, output)
        else:
            export_csv(sheets, output)
        return sheets

    def _resolve_remote_identity(self, mol, decisions: list[FetchDecision]):
        if not decisions or not self.remote_enabled:
            return None
        identity_decisions = [d for d in decisions if d.provider == "PubChemPugRest"]
        if not identity_decisions:
            return None
        decision = identity_decisions[0]
        hits = []
        if decision.job_type == "formula_identity":
            hits = self.pubchem.resolve_formula(
                mol,
                max_cids_per_formula=int(self.config["generation"].get("max_cids_per_formula", 10)),
                allowed_elements=set(self.config["prefilter"].get("allowed_elements", [])) or None,
                max_heavy_atoms=int(self.config["prefilter"].get("max_heavy_atoms", 20)),
                max_molecular_weight=float(self.config["prefilter"].get("max_molecular_weight", 300)),
            )
        elif decision.job_type == "identity":
            hits = self.pubchem.resolve(mol)
        resolved = [h for h in hits if h.identity_status == "resolved"]
        return resolved[0] if resolved else None

    def _classify_pugview(self, chemical, decisions: list[FetchDecision]):
        if not decisions or not any(d.provider == "PubChemPugView" for d in decisions):
            return []
        return self.pugview.classify(chemical)

    def _plan_remote(self, mol, selected) -> list[FetchDecision]:
        if self.max_api_candidates_per_run <= 0:
            return []
        return self.fetch_planner.plan_for_molecule(
            mol,
            selected_properties=selected,
            remote_enabled=self.remote_enabled,
            pubchem_enabled=self.pubchem_enabled,
            pugview_enabled=self.pugview_enabled,
            mode=self.mode,
            allow_formula_search=self._formula_search_enabled(),
        )

    def _planned_api_rows(self, molecules, apply_limit: bool = False) -> list[dict]:
        rows: list[dict] = []
        remote_candidates = 0
        for mol in molecules:
            decisions = self.fetch_planner.plan_for_molecule(
                mol,
                selected_properties=None,
                remote_enabled=self.remote_enabled,
                pubchem_enabled=self.pubchem_enabled,
                pugview_enabled=self.pugview_enabled,
                mode=self.mode,
                allow_formula_search=self._formula_search_enabled(),
            )
            if decisions and apply_limit:
                if remote_candidates >= self.max_api_candidates_per_run:
                    decisions = []
                else:
                    remote_candidates += 1
            rows.extend(self._decision_rows(mol.candidate_id, decisions))
        return rows

    def _decision_rows(self, candidate_id: str, decisions: list[FetchDecision]) -> list[dict]:
        return [
            {
                "candidate_id": candidate_id,
                "provider": decision.provider,
                "job_type": decision.job_type,
                "priority": decision.priority,
                "reason": decision.reason,
            }
            for decision in decisions
        ]

    def _record_planned_api_stats(self, planned_rows: list[dict]) -> None:
        self.stats["planned_api_requests"] = len(planned_rows)
        self.stats["planned_api_candidates"] = len({row["candidate_id"] for row in planned_rows})

    def _formula_search_enabled(self) -> bool:
        return (
            self.mode == "exploration"
            and self.remote_enabled
            and bool(self.config["generation"].get("formula_search_enabled", False))
        )

    def _fetch_properties(self, chemical):
        required = list(self.config["selection"].get("required_properties", [])) + ["vapor_pressure"]
        out = []
        for provider in self.property_providers:
            props = [p for p in required if provider.supports(p)]
            if not props:
                continue
            out.extend(provider.fetch(chemical, props))
        return out

    def _select_properties(self, chemical, candidates):
        required = list(self.config["selection"].get("required_properties", []))
        selected = self.selector.select(candidates, required)
        selected = [s if s.chemical_id else replace(s, chemical_id=chemical.chemical_id) for s in selected]
        for T in self.config["selection"].get("target_vapor_temperatures_K", [298.15, 313.15, 333.15]):
            selected.append(self.selector.select_by_temperature(candidates, "vapor_pressure", float(T)))
        return [s if s.chemical_id else replace(s, chemical_id=chemical.chemical_id) for s in selected]

    def _sheets(self, summary_rows, evidence_rows, coverage_rows, review_rows, rejected_rows, planned_rows):
        cache_stats = self.cache.stats
        self.stats["pubchem_requests"] = cache_stats.get("pubchem_requests", 0)
        self.stats["pugview_requests"] = cache_stats.get("pugview_requests", 0)
        self.stats["cache_hits"] = cache_stats.get("cache_hits", 0)
        self.stats["negative_cache_hits"] = cache_stats.get("negative_cache_hits", 0)
        return {
            "Summary": summary_rows,
            "Evidence": evidence_rows,
            "Coverage": coverage_rows,
            "Review Required": review_rows,
            "Rejected": rejected_rows,
            "Run Stats": self._stats_rows(),
            "Planned API": planned_rows,
            "Candidate Breakdown": candidate_breakdown_rows(summary_rows),
            "Amine Summary": amine_summary_rows(summary_rows),
            "Amine Class Summary": amine_class_summary_rows(summary_rows),
        }

    def _empty_stats(self) -> dict[str, int]:
        return {
            "generated_candidates": 0,
            "normalized_candidates": 0,
            "deduplicated_candidates": 0,
            "prefilter_passed": 0,
            "prefilter_rejected": 0,
            "pubchem_requests": 0,
            "pugview_requests": 0,
            "cache_hits": 0,
            "negative_cache_hits": 0,
            "planned_api_requests": 0,
            "planned_api_candidates": 0,
            "remote_enrichment_candidates": 0,
            "final_summary_rows": 0,
        }

    def _stats_rows(self) -> list[dict]:
        return [{"metric": key, "value": value} for key, value in self.stats.items()]

    def _evidence_rows(self, chemical, candidates, selected, classifications, kinetics):
        rows = []
        for c in candidates:
            rows.append({
                "chemical_id": chemical.chemical_id,
                "preferred_name": chemical.preferred_name,
                "record_type": "property_candidate",
                "property_name": c.property_name,
                "value_num": c.value_num,
                "value_text": c.value_text,
                "unit": c.unit,
                "temperature_K": c.temperature_K,
                "pressure_Pa": c.pressure_Pa,
                "source": c.source,
                "method": c.method,
                "is_estimated": c.is_estimated,
                "quality_hint": c.quality_hint,
                "reference": c.reference,
                "selected": False,
                "selection_reason": "",
            })
        for s in selected:
            rows.append({
                "chemical_id": chemical.chemical_id,
                "preferred_name": chemical.preferred_name,
                "record_type": "selected_property",
                "property_name": s.property_name,
                "value_num": s.value_num,
                "value_text": s.value_text,
                "unit": s.unit,
                "temperature_K": None,
                "pressure_Pa": None,
                "source": s.selected_source,
                "method": None,
                "is_estimated": None,
                "quality_hint": s.quality_rank,
                "reference": None,
                "selected": True,
                "selection_reason": s.selection_reason,
            })
        for r in classifications:
            rows.append({
                "chemical_id": chemical.chemical_id,
                "preferred_name": chemical.preferred_name,
                "record_type": "classification",
                "property_name": r.classification_type,
                "value_num": None,
                "value_text": r.value,
                "unit": None,
                "temperature_K": None,
                "pressure_Pa": None,
                "source": r.source,
                "method": r.basis,
                "is_estimated": None,
                "quality_hint": r.confidence,
                "reference": None,
                "selected": True,
                "selection_reason": r.basis,
            })
        for k in kinetics:
            rows.append({
                "chemical_id": chemical.chemical_id,
                "preferred_name": chemical.preferred_name,
                "record_type": "kinetics_probe",
                "property_name": f"{k.target_species}:{k.target_species_state or ''}",
                "value_num": None,
                "value_text": k.availability,
                "unit": None,
                "temperature_K": None,
                "pressure_Pa": None,
                "source": "; ".join(k.sources),
                "method": k.phase,
                "is_estimated": None,
                "quality_hint": None,
                "reference": None,
                "selected": True,
                "selection_reason": k.comment,
            })
        return rows

    def _review_row(self, summary):
        issues = []
        if summary.get("identity_status") != "resolved":
            issues.append("identity_not_resolved")
        if summary.get("data_quality") in {"Conflict", "Partial", "Missing", "D"}:
            issues.append(f"data_quality={summary.get('data_quality')}")
        if summary.get("pfas_flag") in {"yes", "possible", "unknown"}:
            issues.append(f"pfas_flag={summary.get('pfas_flag')}")
        if summary.get("kinetics_coverage") == "partial":
            issues.append("kinetics_partial")
        return {
            "preferred_name": summary.get("preferred_name"),
            "cas": summary.get("cas"),
            "issue_type": "; ".join(issues),
            "recommended_action": "Review identity, property evidence, PFAS basis, and shortlisted reaction data before engineering release.",
        }

    def _coverage_rows(self, summary_rows):
        if not summary_rows:
            return []
        props = [
            "tm_C", "tb_C", "tc_C", "pc_MPa", "pvap_25C_kPa", "pvap_40C_kPa", "pvap_60C_kPa",
            "gwp100_ar6", "pfas_flag", "persistence_screen", "reactive_groups", "kinetics_coverage",
        ]
        rows = []
        total = len(summary_rows)
        for p in props:
            present = sum(1 for r in summary_rows if r.get(p) not in (None, "", "unknown", "not_checked"))
            rows.append({
                "property_name": p,
                "total_molecules": total,
                "present_or_informative_count": present,
                "missing_or_unknown_count": total - present,
                "coverage_rate": present / total if total else 0,
            })
        return rows


def run_pipeline(
    config_path: str | Path | None = None,
    input_csv: str | Path | None = None,
    output_xlsx: str | Path | None = None,
    mode: str | None = None,
    remote: bool | None = None,
    max_candidates: int | None = None,
    dry_run: bool = False,
):
    cfg = load_config(config_path)
    _apply_runtime_overrides(cfg, mode=mode, remote=remote, max_candidates=max_candidates)
    pipe = ScreeningPipeline(cfg)
    if output_xlsx:
        return pipe.run_to_output(output_xlsx, input_csv=input_csv, dry_run=dry_run)
    return pipe.run(input_csv=input_csv, dry_run=dry_run)


def _apply_runtime_overrides(cfg: dict, mode: str | None = None, remote: bool | None = None, max_candidates: int | None = None) -> None:
    if mode is not None:
        cfg["mode"] = _normalize_mode(mode)
    if remote is not None:
        cfg["providers"]["pubchem_enabled"] = bool(remote)
        cfg["providers"]["pugview_enabled"] = bool(remote)
    if max_candidates is not None:
        cfg["generation"]["max_total_candidates"] = int(max_candidates)


def _normalize_mode(mode: str | None) -> str:
    value = (mode or "enrichment").strip().lower()
    if value not in MODES:
        raise ValueError(f"Unsupported mode: {mode!r}. Expected one of {sorted(MODES)}")
    return value
