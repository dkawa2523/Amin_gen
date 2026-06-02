from __future__ import annotations

from collections import defaultdict

from gas_screening_mvp.domain.models import PropertyCandidate, SelectedProperty
from gas_screening_mvp.selection.unit_normalizer import normalize_property


DEFAULT_PROVIDER_PRIORITY = {
    "CuratedCsv": 5,
    "GWP_local": 5,
    "CoolProp": 10,
    "chemicals": 20,
    "thermo.VaporPressure": 30,
    "LocalAmineEstimate": 70,
    "PubChemPugView": 80,
}

DEFAULT_TARGET_UNITS = {
    "normal_boiling_point": "K",
    "normal_melting_point": "K",
    "critical_temperature": "K",
    "critical_pressure": "Pa",
    "vapor_pressure": "Pa",
    "gwp100_ar6": "kg_CO2e_per_kg",
    "gwp100_ar5": "kg_CO2e_per_kg",
    "gwp100_ar4": "kg_CO2e_per_kg",
}

CONFLICT_THRESHOLDS = {
    "normal_boiling_point": 3.0,  # K
    "normal_melting_point": 5.0,
    "critical_temperature": 5.0,
    "critical_pressure": 0.05,  # relative fraction if values large
    "vapor_pressure": 0.30,     # relative fraction
}


class PropertySelector:
    def __init__(self, provider_priority: dict[str, int] | None = None, target_units: dict[str, str] | None = None):
        self.provider_priority = dict(DEFAULT_PROVIDER_PRIORITY)
        if provider_priority:
            self.provider_priority.update(provider_priority)
        self.target_units = dict(DEFAULT_TARGET_UNITS)
        if target_units:
            self.target_units.update(target_units)

    def select(self, candidates: list[PropertyCandidate], required_properties: list[str]) -> list[SelectedProperty]:
        grouped: dict[str, list[PropertyCandidate]] = defaultdict(list)
        for c in candidates:
            grouped[c.property_name].append(normalize_property(c, self.target_units.get(c.property_name)))

        out: list[SelectedProperty] = []
        for pname in required_properties:
            items = grouped.get(pname, [])
            if not items:
                out.append(SelectedProperty(
                    chemical_id=candidates[0].chemical_id if candidates else "",
                    property_name=pname,
                    value_num=None,
                    value_text=None,
                    unit=self.target_units.get(pname),
                    status="missing",
                    quality_rank="Missing",
                    selected_source=None,
                    selection_reason="No candidate values available",
                ))
                continue
            out.append(self._select_one(pname, items))
        return out

    def select_by_temperature(self, candidates: list[PropertyCandidate], property_name: str, target_T_K: float) -> SelectedProperty:
        items = [normalize_property(c, self.target_units.get(property_name)) for c in candidates if c.property_name == property_name]
        same_T = [c for c in items if c.temperature_K is not None and abs(c.temperature_K - target_T_K) < 0.2]
        if not same_T:
            return SelectedProperty(
                chemical_id=candidates[0].chemical_id if candidates else "",
                property_name=f"{property_name}_{round(target_T_K,2)}K",
                value_num=None,
                value_text=None,
                unit=self.target_units.get(property_name),
                status="missing",
                quality_rank="Missing",
                selected_source=None,
                selection_reason="No vapor pressure candidate at target temperature",
            )
        ranked = sorted([c for c in same_T if c.value_num is not None or c.value_text is not None], key=self._rank_key)
        best_candidate = ranked[0] if ranked else None
        selected = self._select_one(property_name, same_T)
        status = selected.status
        quality_rank = selected.quality_rank
        selection_reason = selected.selection_reason
        if selected.status == "selected" and best_candidate is not None:
            outside = _outside_temperature_range(
                target_T_K,
                best_candidate.valid_temperature_min_K,
                best_candidate.valid_temperature_max_K,
            )
            if outside:
                status = "outside_range"
                quality_rank = "D"
                selection_reason = (
                    "Selected highest-priority source, but target temperature is outside "
                    f"the candidate valid range: {outside}"
                )
        return SelectedProperty(
            chemical_id=selected.chemical_id,
            property_name=f"{property_name}_{round(target_T_K,2)}K",
            value_num=selected.value_num,
            value_text=selected.value_text,
            unit=selected.unit,
            status=status,
            quality_rank=quality_rank,
            selected_source=selected.selected_source,
            selection_reason=selection_reason,
        )

    def _select_one(self, pname: str, items: list[PropertyCandidate]) -> SelectedProperty:
        items = [c for c in items if c.value_num is not None or c.value_text is not None]
        if not items:
            return SelectedProperty("", pname, None, None, self.target_units.get(pname), "missing", "Missing", None, "Only empty candidate records")

        sorted_items = sorted(items, key=self._rank_key)
        best = sorted_items[0]
        conflict = self._detect_conflict(pname, items)
        if conflict:
            return SelectedProperty(
                best.chemical_id,
                pname,
                best.value_num,
                best.value_text,
                best.unit,
                "conflict",
                "Conflict",
                best.source,
                "Selected highest-priority source but conflict detected across candidate values",
            )
        quality = self._quality(best)
        return SelectedProperty(
            best.chemical_id,
            pname,
            best.value_num,
            best.value_text,
            best.unit,
            "selected",
            quality,
            best.source,
            f"Selected by provider priority and quality: {best.source}/{best.method or ''}",
        )

    def _rank_key(self, c: PropertyCandidate):
        source_rank = self.provider_priority.get(c.source, 50)
        est_penalty = 30 if c.is_estimated else 0
        quality_penalty = {"A": 0, "B": 5, "C": 15, "D": 30}.get(c.quality_hint or "C", 15)
        return (source_rank + est_penalty + quality_penalty, c.retrieved_at)

    def _quality(self, c: PropertyCandidate) -> str:
        if c.quality_hint in {"A", "B", "C", "D"}:
            return c.quality_hint
        if c.is_estimated:
            return "C"
        return "B"

    def _detect_conflict(self, pname: str, items: list[PropertyCandidate]) -> bool:
        comparable_items = [c for c in items if not c.is_estimated]
        if len(comparable_items) == 1:
            return False
        if not comparable_items:
            comparable_items = items
        vals = [c.value_num for c in comparable_items if c.value_num is not None]
        if len(vals) < 2:
            return False
        lo, hi = min(vals), max(vals)
        threshold = CONFLICT_THRESHOLDS.get(pname)
        if threshold is None:
            return False
        if pname in {"vapor_pressure", "critical_pressure"}:
            if hi == 0:
                return False
            return (hi - lo) / hi > threshold
        return (hi - lo) > threshold


def _outside_temperature_range(T_K: float, valid_min_K: float | None, valid_max_K: float | None) -> str | None:
    if valid_min_K is not None and T_K < valid_min_K:
        return f"{T_K:g} K < {valid_min_K:g} K"
    if valid_max_K is not None and T_K > valid_max_K:
        return f"{T_K:g} K > {valid_max_K:g} K"
    return None
