from __future__ import annotations

from typing import Iterable

from gas_screening_mvp.domain.models import ChemicalIdentity, PropertyCandidate


class ThermoVaporPressureProvider:
    """Optional vapor pressure provider using `thermo.VaporPressure`.

    It returns vapor pressure candidates at configured target temperatures.
    Extrapolation is not explicitly requested; invalid results are ignored.
    """

    name = "thermo.VaporPressure"

    def __init__(self, target_temperatures_K: list[float] | None = None):
        self.target_temperatures_K = target_temperatures_K or [298.15, 313.15, 333.15]

    def supports(self, property_name: str) -> bool:
        return property_name == "vapor_pressure"

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        if "vapor_pressure" not in set(property_names) or not chemical.cas:
            return []
        try:
            from thermo import VaporPressure
        except Exception:
            return []
        try:
            vp = VaporPressure(CASRN=chemical.cas)
        except Exception:
            return []
        out: list[PropertyCandidate] = []
        for T in self.target_temperatures_K:
            try:
                value = vp.T_dependent_property(T)
            except Exception:
                value = None
            if value is None:
                continue
            out.append(
                PropertyCandidate(
                    chemical_id=chemical.chemical_id,
                    property_name="vapor_pressure",
                    value_num=float(value),
                    unit="Pa",
                    temperature_K=T,
                    source=self.name,
                    method=getattr(vp, "method", None),
                    valid_temperature_min_K=getattr(vp, "Tmin", None),
                    valid_temperature_max_K=getattr(vp, "Tmax", None),
                    is_estimated=False,
                    quality_hint="B",
                )
            )
        return out
