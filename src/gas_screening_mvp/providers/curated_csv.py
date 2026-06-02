from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from gas_screening_mvp.domain.models import ChemicalIdentity, PropertyCandidate


class CuratedCsvPropertyProvider:
    """Loads property candidates from a CSV table.

    Expected columns:
      key_type,key_value,property_name,value_num,value_text,unit,temperature_K,
      pressure_Pa,source,method,is_estimated,quality_hint,reference,
      valid_temperature_min_K,valid_temperature_max_K

    key_type can be inchikey, cas, name, pubchem_cid.
    """

    name = "CuratedCsv"

    def __init__(self, csv_path: str | Path | None):
        self.csv_path = Path(csv_path) if csv_path else None
        self.rows: list[dict] = []
        if self.csv_path and self.csv_path.exists():
            with self.csv_path.open(newline="", encoding="utf-8-sig") as f:
                self.rows = list(csv.DictReader(f))

    def supports(self, property_name: str) -> bool:
        return True

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        wanted = set(property_names)
        out: list[PropertyCandidate] = []
        keys = {
            "inchikey": chemical.inchikey,
            "cas": chemical.cas,
            "name": (chemical.preferred_name or "").lower() if chemical.preferred_name else None,
            "pubchem_cid": str(chemical.pubchem_cid) if chemical.pubchem_cid is not None else None,
        }
        for row in self.rows:
            pname = row.get("property_name")
            if wanted and pname not in wanted:
                continue
            kt = (row.get("key_type") or "").strip().lower()
            kv = (row.get("key_value") or "").strip()
            if kt == "name":
                match = keys.get(kt) == kv.lower()
            else:
                match = keys.get(kt) == kv
            if not match:
                continue
            out.append(
                PropertyCandidate(
                    chemical_id=chemical.chemical_id,
                    property_name=pname,
                    value_num=_float_or_none(row.get("value_num")),
                    value_text=row.get("value_text") or None,
                    unit=row.get("unit") or None,
                    temperature_K=_float_or_none(row.get("temperature_K")),
                    pressure_Pa=_float_or_none(row.get("pressure_Pa")),
                    source=row.get("source") or self.name,
                    method=row.get("method") or "curated_csv",
                    is_estimated=_bool(row.get("is_estimated")),
                    quality_hint=row.get("quality_hint") or None,
                    reference=row.get("reference") or None,
                    valid_temperature_min_K=_float_or_none(row.get("valid_temperature_min_K")),
                    valid_temperature_max_K=_float_or_none(row.get("valid_temperature_max_K")),
                )
            )
        return out


def _float_or_none(v):
    try:
        if v is None or v == "":
            return None
        return float(v)
    except Exception:
        return None


def _bool(v) -> bool:
    return str(v).strip().lower() in {"1", "true", "yes", "y"}
