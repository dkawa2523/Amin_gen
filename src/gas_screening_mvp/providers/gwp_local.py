from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from gas_screening_mvp.domain.models import ChemicalIdentity, PropertyCandidate


class GwpLocalProvider:
    """Loads local GWP records.

    Expected columns include any of: inchikey, cas, name, formula,
    gwp100_ar6, gwp100_ar5, gwp100_ar4, source_version, reference.
    """

    name = "GWP_local"

    def __init__(self, csv_path: str | Path | None):
        self.rows = []
        self.csv_path = Path(csv_path) if csv_path else None
        if self.csv_path and self.csv_path.exists():
            with self.csv_path.open(newline="", encoding="utf-8-sig") as f:
                self.rows = list(csv.DictReader(f))

    def supports(self, property_name: str) -> bool:
        return property_name in {"gwp100_ar6", "gwp100_ar5", "gwp100_ar4"}

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        wanted = set(property_names)
        out: list[PropertyCandidate] = []
        for row in self.rows:
            if not self._matches(row, chemical):
                continue
            for pname in wanted:
                val = row.get(pname)
                if val not in (None, ""):
                    out.append(PropertyCandidate(
                        chemical_id=chemical.chemical_id,
                        property_name=pname,
                        value_num=_float_or_none(val),
                        unit="kg_CO2e_per_kg",
                        source=self.name,
                        method="local_table",
                        source_version=row.get("source_version") or None,
                        reference=row.get("reference") or None,
                        is_estimated=False,
                        quality_hint="A",
                    ))
        return out

    def _matches(self, row: dict, chemical: ChemicalIdentity) -> bool:
        if row.get("inchikey") and chemical.inchikey and row["inchikey"] == chemical.inchikey:
            return True
        if row.get("cas") and chemical.cas and row["cas"] == chemical.cas:
            return True
        if row.get("formula") and chemical.formula and row["formula"].replace(" ", "") == chemical.formula.replace(" ", ""):
            return True
        if row.get("name") and chemical.preferred_name and row["name"].lower() == chemical.preferred_name.lower():
            return True
        return False


def _float_or_none(v):
    try:
        return float(v)
    except Exception:
        return None
