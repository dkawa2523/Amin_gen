from __future__ import annotations

from typing import Iterable

from gas_screening_mvp.domain.models import ChemicalIdentity, PropertyCandidate


class ChemicalsLocalProvider:
    """Optional provider using the `chemicals` Python package.

    If the dependency is not installed or CAS is missing, it returns no data.
    """

    name = "chemicals"

    PROPERTY_TO_FUNC = {
        "normal_boiling_point": ("chemicals.phase_change", "Tb", "K"),
        "normal_melting_point": ("chemicals.phase_change", "Tm", "K"),
        "critical_temperature": ("chemicals.critical", "Tc", "K"),
        "critical_pressure": ("chemicals.critical", "Pc", "Pa"),
    }

    def supports(self, property_name: str) -> bool:
        return property_name in self.PROPERTY_TO_FUNC

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        if not chemical.cas:
            return []
        out: list[PropertyCandidate] = []
        for pname in property_names:
            if pname not in self.PROPERTY_TO_FUNC:
                continue
            mod_name, func_name, unit = self.PROPERTY_TO_FUNC[pname]
            try:
                mod = __import__(mod_name, fromlist=[func_name])
                func = getattr(mod, func_name)
                value = func(CASRN=chemical.cas)
            except Exception:
                value = None
            if value is not None:
                out.append(
                    PropertyCandidate(
                        chemical_id=chemical.chemical_id,
                        property_name=pname,
                        value_num=float(value),
                        unit=unit,
                        source=self.name,
                        method=func_name,
                        is_estimated=False,
                        quality_hint="A",
                    )
                )
        return out
