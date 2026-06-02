from __future__ import annotations

from typing import Iterable

from gas_screening_mvp.domain.models import ChemicalIdentity, PropertyCandidate


COOLPROP_NAME_MAP = {
    "nitrogen": "Nitrogen",
    "oxygen": "Oxygen",
    "argon": "Argon",
    "helium": "Helium",
    "hydrogen": "Hydrogen",
    "carbon dioxide": "CarbonDioxide",
    "carbon monoxide": "CarbonMonoxide",
    "ammonia": "Ammonia",
    "methane": "Methane",
    "sulfur hexafluoride": "SulfurHexafluoride",
    # Extend in deployment with validated CoolProp names.
}


class CoolPropLocalProvider:
    name = "CoolProp"

    def __init__(self, target_temperatures_K: list[float] | None = None):
        self.target_temperatures_K = target_temperatures_K or [298.15, 313.15, 333.15]

    def supports(self, property_name: str) -> bool:
        return property_name in {"critical_temperature", "critical_pressure", "vapor_pressure"}

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        wanted = set(property_names)
        fluid = self._fluid_name(chemical)
        if not fluid:
            return []
        try:
            from CoolProp.CoolProp import PropsSI
        except Exception:
            return []
        out: list[PropertyCandidate] = []
        if "critical_temperature" in wanted:
            try:
                out.append(PropertyCandidate(chemical.chemical_id, "critical_temperature", float(PropsSI("Tcrit", fluid)), None, "K", source=self.name, method="PropsSI", quality_hint="A"))
            except Exception:
                pass
        if "critical_pressure" in wanted:
            try:
                out.append(PropertyCandidate(chemical.chemical_id, "critical_pressure", float(PropsSI("pcrit", fluid)), None, "Pa", source=self.name, method="PropsSI", quality_hint="A"))
            except Exception:
                pass
        if "vapor_pressure" in wanted:
            for T in self.target_temperatures_K:
                try:
                    p = PropsSI("P", "T", T, "Q", 0, fluid)
                except Exception:
                    p = None
                if p is not None:
                    out.append(PropertyCandidate(chemical.chemical_id, "vapor_pressure", float(p), None, "Pa", temperature_K=T, source=self.name, method="PropsSI_Q0", quality_hint="A"))
        return out

    def _fluid_name(self, chemical: ChemicalIdentity) -> str | None:
        names = [chemical.preferred_name]
        if chemical.formula:
            names.append(chemical.formula)
        for n in names:
            if not n:
                continue
            key = n.lower().replace("_", " ").strip()
            if key in COOLPROP_NAME_MAP:
                return COOLPROP_NAME_MAP[key]
        return None
