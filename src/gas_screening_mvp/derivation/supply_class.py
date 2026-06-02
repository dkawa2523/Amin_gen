from __future__ import annotations


def derive_supply_class(phase_25C: str, pvap_25_kPa: float | None, pvap_60_kPa: float | None, thresholds: dict | None = None) -> str:
    th = thresholds or {}
    high = float(th.get("high_vapor_pressure_kPa", 100.0))
    bubbler = float(th.get("bubbler_kPa", 10.0))
    heated = float(th.get("heated_source_kPa", 1.0))

    if phase_25C in {"gas", "gas_or_supercritical"}:
        return "compressed_or_liquefied_gas"
    if pvap_25_kPa is not None:
        if pvap_25_kPa >= high:
            return "high_vapor_pressure_liquid"
        if pvap_25_kPa >= bubbler:
            return "bubbler_or_direct_liquid_source"
        if pvap_25_kPa >= heated:
            return "heated_source_likely"
    if pvap_60_kPa is not None and pvap_60_kPa >= heated:
        return "heated_source_required"
    if phase_25C == "solid":
        return "solid_or_sublimation_source_review"
    if phase_25C == "unknown":
        return "unknown_review_required"
    return "low_volatility_or_special_vaporizer"
