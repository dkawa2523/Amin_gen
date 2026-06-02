from __future__ import annotations

from dataclasses import dataclass

from gas_screening_mvp.domain.models import NormalizedMolecule, SelectedProperty
from gas_screening_mvp.prefilter.score import api_priority_score
from gas_screening_mvp.storage.cache import SqliteApiCache


@dataclass(frozen=True)
class FetchDecision:
    candidate_id: str
    provider: str
    job_type: str
    priority: int
    reason: str


class FetchPlanner:
    """Decides which remote providers should be called.

    The MVP uses this for gating remote enrichment. The actual provider calls are
    orchestrated by the pipeline.
    """

    def __init__(self, cache: SqliteApiCache, thresholds: dict | None = None):
        self.cache = cache
        self.thresholds = thresholds or {}
        self.identity_threshold = float(self.thresholds.get("min_identity_api_score", 0.45))
        self.enrichment_threshold = float(self.thresholds.get("min_enrichment_api_score", 0.85))

    def plan_for_molecule(
        self,
        mol: NormalizedMolecule,
        selected_properties: list[SelectedProperty] | None = None,
        remote_enabled: bool = True,
        pubchem_enabled: bool | None = None,
        pugview_enabled: bool | None = None,
        mode: str = "enrichment",
        allow_formula_search: bool = False,
        shortlist: bool = False,
    ) -> list[FetchDecision]:
        if not remote_enabled:
            return []
        pubchem_allowed = remote_enabled if pubchem_enabled is None else bool(pubchem_enabled)
        pugview_allowed = remote_enabled if pugview_enabled is None else bool(pugview_enabled)
        score = api_priority_score(mol)
        decisions: list[FetchDecision] = []
        has_missing = _has_missing_or_weak_selection(selected_properties)
        if not has_missing:
            if pugview_allowed and (shortlist or score >= self.enrichment_threshold):
                decisions.append(FetchDecision(mol.candidate_id, "PubChemPugView", "ghs", int(score * 100), "shortlist_or_score_above_enrichment_threshold"))
            return decisions

        if score < self.identity_threshold:
            return decisions

        if pubchem_allowed and mol.metadata.get("pubchem_cid"):
            decisions.append(FetchDecision(mol.candidate_id, "PubChemPugRest", "identity", int(score * 100), "pubchem_cid_from_generation"))
        elif pubchem_allowed and mol.standard_inchikey and not self.cache.has_negative("PubChemPugRest", "inchikey", mol.standard_inchikey):
            decisions.append(FetchDecision(mol.candidate_id, "PubChemPugRest", "identity", int(score * 100), "missing_data_and_score_above_identity_threshold"))
        elif pubchem_allowed and mol.canonical_smiles and not self.cache.has_negative("PubChemPugRest", "smiles", mol.canonical_smiles):
            decisions.append(FetchDecision(mol.candidate_id, "PubChemPugRest", "identity", int(score * 100), "missing_data_and_score_above_identity_threshold"))
        elif (
            pubchem_allowed
            and mode == "exploration"
            and allow_formula_search
            and mol.formula
            and not mol.canonical_smiles
            and not self.cache.has_negative("PubChemPugRest", "formula", mol.formula)
        ):
            decisions.append(FetchDecision(mol.candidate_id, "PubChemPugRest", "formula_identity", int(score * 100), "formula_only_exploration_search"))

        if pugview_allowed and (shortlist or score >= self.enrichment_threshold):
            decisions.append(FetchDecision(mol.candidate_id, "PubChemPugView", "ghs", int(score * 100), "shortlist_or_score_above_enrichment_threshold"))
        return decisions


def _has_missing_or_weak_selection(selected_properties: list[SelectedProperty] | None) -> bool:
    if selected_properties is None:
        return True
    if not selected_properties:
        return True
    weak_status = {"missing", "conflict", "outside_range", "manual_review_required"}
    return any(prop.status in weak_status or prop.quality_rank in {"Missing", "Conflict", "D"} for prop in selected_properties)
