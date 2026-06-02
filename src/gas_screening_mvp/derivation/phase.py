from __future__ import annotations


def pvap_status(T_use_K: float, Tc_K: float | None, selected_status: str, valid_min_K: float | None = None, valid_max_K: float | None = None) -> str:
    if Tc_K is not None and T_use_K >= Tc_K:
        return "not_applicable_above_Tc"
    if selected_status == "missing":
        return "no_correlation"
    if valid_min_K is not None and T_use_K < valid_min_K:
        return "outside_correlation_range"
    if valid_max_K is not None and T_use_K > valid_max_K:
        return "outside_correlation_range"
    if selected_status == "selected":
        return "value_available"
    return selected_status


def phase_at_25C(Tm_K: float | None, Tb_K: float | None, Tc_K: float | None, T_K: float = 298.15) -> str:
    if Tc_K is not None and T_K >= Tc_K:
        return "gas_or_supercritical"
    if Tm_K is None or Tb_K is None:
        return "unknown"
    if T_K < Tm_K:
        return "solid"
    if Tm_K <= T_K < Tb_K:
        return "liquid"
    return "gas"
