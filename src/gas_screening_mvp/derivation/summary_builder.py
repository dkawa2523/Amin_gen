from __future__ import annotations

from collections import defaultdict

from gas_screening_mvp.domain.models import ChemicalIdentity, SelectedProperty, ClassificationRecord, ReactionProbeSummary
from gas_screening_mvp.derivation.phase import phase_at_25C, pvap_status
from gas_screening_mvp.derivation.supply_class import derive_supply_class


def selected_map(selected: list[SelectedProperty]) -> dict[str, SelectedProperty]:
    return {s.property_name: s for s in selected}


def _val(sel: dict[str, SelectedProperty], name: str):
    s = sel.get(name)
    return s.value_num if s else None


def _to_C(K: float | None) -> float | None:
    return None if K is None else K - 273.15


def _pa_to_kpa(Pa: float | None) -> float | None:
    return None if Pa is None else Pa / 1000.0


def _pa_to_mpa(Pa: float | None) -> float | None:
    return None if Pa is None else Pa / 1_000_000.0


def _join(values) -> str:
    vals = [v for v in values if v]
    return "; ".join(sorted(set(vals))) if vals else ""


def build_summary_row(
    chemical: ChemicalIdentity,
    selected: list[SelectedProperty],
    classifications: list[ClassificationRecord],
    kinetics: list[ReactionProbeSummary],
    thresholds: dict | None = None,
) -> dict[str, object]:
    sel = selected_map(selected)
    Tm = _val(sel, "normal_melting_point")
    Tb = _val(sel, "normal_boiling_point")
    Tc = _val(sel, "critical_temperature")
    Pc = _val(sel, "critical_pressure")
    pvap25 = _val(sel, "vapor_pressure_298.15K")
    pvap40 = _val(sel, "vapor_pressure_313.15K")
    pvap60 = _val(sel, "vapor_pressure_333.15K")
    pvap25_sel = sel.get("vapor_pressure_298.15K")
    pvap40_sel = sel.get("vapor_pressure_313.15K")
    pvap60_sel = sel.get("vapor_pressure_333.15K")
    gwp100_ar6 = sel.get("gwp100_ar6")

    phase = phase_at_25C(Tm, Tb, Tc)
    supply = derive_supply_class(phase, _pa_to_kpa(pvap25), _pa_to_kpa(pvap60), thresholds)

    cls = defaultdict(list)
    for rec in classifications:
        cls[rec.classification_type].append(rec.value)

    pfas_flag = cls.get("pfas_flag", ["unknown"])[0]
    pfas_basis = cls.get("pfas_basis", [""])[0]
    persistence_screen = cls.get("persistence_screen", ["unknown"])[0]
    persistence_basis = cls.get("persistence_basis", [""])[0]

    avail_order = {"available": 3, "partial": 2, "none": 1, "not_checked": 0}
    best_cov = "not_checked"
    sources = []
    for k in kinetics:
        if avail_order.get(k.availability, 0) > avail_order.get(best_cov, 0):
            best_cov = k.availability
        sources.extend(k.sources)

    quality_values = [s.quality_rank for s in selected]
    if "Conflict" in quality_values:
        data_quality = "Conflict"
    elif "Missing" in quality_values:
        data_quality = "Partial"
    elif quality_values:
        rank_order = {"A": 1, "B": 2, "C": 3, "D": 4, "N/A": 5}
        data_quality = max(quality_values, key=lambda x: rank_order.get(x, 6))
    else:
        data_quality = "Missing"

    review_required = bool(
        chemical.identity_status != "resolved"
        or data_quality in {"Conflict", "Partial", "Missing", "D"}
        or pfas_flag in {"yes", "possible", "unknown"}
        or best_cov == "partial"
    )

    return {
        "input_name": chemical.preferred_name,
        "preferred_name": chemical.preferred_name,
        "cas": chemical.cas,
        "pubchem_cid": chemical.pubchem_cid,
        "formula": chemical.formula,
        "molecular_weight": chemical.molecular_weight,
        "canonical_smiles": chemical.canonical_smiles,
        "inchikey": chemical.inchikey,
        "identity_status": chemical.identity_status,
        "tm_C": _to_C(Tm),
        "tb_C": _to_C(Tb),
        "tc_C": _to_C(Tc),
        "pc_MPa": _pa_to_mpa(Pc),
        "pvap_25C_kPa": _pa_to_kpa(pvap25),
        "pvap_25C_status": pvap_status(298.15, Tc, pvap25_sel.status if pvap25_sel else "missing"),
        "pvap_40C_kPa": _pa_to_kpa(pvap40),
        "pvap_40C_status": pvap_status(313.15, Tc, pvap40_sel.status if pvap40_sel else "missing"),
        "pvap_60C_kPa": _pa_to_kpa(pvap60),
        "pvap_60C_status": pvap_status(333.15, Tc, pvap60_sel.status if pvap60_sel else "missing"),
        "phase_25C_1atm": phase,
        "supply_class": supply,
        "gwp100_ar6": _val(sel, "gwp100_ar6"),
        "gwp100_ar6_status": gwp100_ar6.status if gwp100_ar6 else "missing",
        "gwp100_ar6_source": gwp100_ar6.selected_source if gwp100_ar6 else None,
        "pfas_flag": pfas_flag,
        "pfas_basis": pfas_basis,
        "pfas_list_hits": _join(cls.get("pfas_list_hit", [])),
        "persistence_screen": persistence_screen,
        "persistence_basis": persistence_basis,
        "reactive_groups": _join(cls.get("reactive_group", [])),
        "reactivity_flags": _join(cls.get("reactivity_flag", [])),
        "ghs_physical_h_codes": _join(cls.get("ghs_physical_h_code", [])),
        "kinetics_coverage": best_cov,
        "kinetics_sources": _join(sources),
        "data_quality": data_quality,
        "review_required": review_required,
    }
