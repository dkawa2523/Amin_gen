from __future__ import annotations

import csv
from pathlib import Path

from gas_screening_mvp.domain.models import ChemicalIdentity, ClassificationRecord

try:
    from rdkit import Chem
except Exception:  # pragma: no cover
    Chem = None


PHYSICAL_H_CODES = {
    "H220", "H221", "H224", "H225", "H226", "H230", "H231",
    "H240", "H241", "H250", "H251", "H252", "H260", "H261",
    "H270", "H271", "H272", "H290",
}


class ReactivityClassifier:
    """CAMEO/local H-code plus RDKit functional group screening."""

    name = "ReactivityClassifier"

    def __init__(self, local_csv: str | Path | None = None):
        self.local_rows: list[dict] = []
        if local_csv and Path(local_csv).exists():
            with Path(local_csv).open(newline="", encoding="utf-8-sig") as f:
                self.local_rows = list(csv.DictReader(f))

    def classify(self, chemical: ChemicalIdentity) -> list[ClassificationRecord]:
        records: list[ClassificationRecord] = []
        records.extend(self._classify_from_local(chemical))
        records.extend(self._classify_from_structure(chemical))
        # Deduplicate by (type,value)
        dedup = {}
        for r in records:
            dedup[(r.classification_type, r.value)] = r
        return list(dedup.values())

    def _classify_from_local(self, chemical: ChemicalIdentity) -> list[ClassificationRecord]:
        out: list[ClassificationRecord] = []
        for row in self.local_rows:
            if not self._matches(row, chemical):
                continue
            for group in _split(row.get("reactive_groups")):
                out.append(ClassificationRecord(chemical.chemical_id, "reactive_group", group, "local curated reactivity table", self.name, confidence="B"))
            for flag in _split(row.get("reactivity_flags")):
                out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", flag, "local curated reactivity table", self.name, confidence="B"))
            for h in _split(row.get("ghs_physical_h_codes")):
                if h in PHYSICAL_H_CODES:
                    out.append(ClassificationRecord(chemical.chemical_id, "ghs_physical_h_code", h, "local curated GHS physical hazard", self.name, confidence="B"))
        return out

    def _classify_from_structure(self, chemical: ChemicalIdentity) -> list[ClassificationRecord]:
        smiles = chemical.canonical_smiles or chemical.isomeric_smiles
        out: list[ClassificationRecord] = []
        if not smiles or Chem is None:
            return out
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return out

        elems = {atom.GetSymbol() for atom in mol.GetAtoms()}
        name = (chemical.preferred_name or "").lower()
        formula = chemical.formula or ""

        # Amines / basic N compounds: N with at least one C neighbor and not only inorganic N.
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7 and any(n.GetAtomicNum() == 6 for n in atom.GetNeighbors()):
                out.append(ClassificationRecord(chemical.chemical_id, "reactive_group", "Amines / basic nitrogen compounds", "RDKit N-C motif", self.name, confidence="C"))
                out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", "acid_base_reactive", "amine/basic nitrogen motif", self.name, confidence="C"))
                break

        if "Si" in elems and ("Cl" in elems or "F" in elems):
            out.append(ClassificationRecord(chemical.chemical_id, "reactive_group", "Silicon halides / halosilanes", "Si-halogen motif", self.name, confidence="C"))
            out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", "water_reactive", "Si-halogen motif", self.name, confidence="C"))

        if "B" in elems and ("Cl" in elems or "F" in elems):
            out.append(ClassificationRecord(chemical.chemical_id, "reactive_group", "Boron halides", "B-halogen motif", self.name, confidence="C"))
            out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", "Lewis_acid", "B-halogen motif", self.name, confidence="C"))

        if "F2" in formula or name == "fluorine" or formula == "F2":
            out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", "strong_oxidizer", "fluorine screening rule", self.name, confidence="C"))
        if name in {"oxygen", "nitrous oxide", "nitrogen trifluoride"}:
            out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", "oxidizer", "known process gas rule", self.name, confidence="C"))
        if name in {"silane", "disilane", "phosphine", "arsine", "hydrogen"}:
            out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", "reducing", "known process gas rule", self.name, confidence="C"))
        if name in {"silane", "disilane", "phosphine", "arsine"}:
            out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", "pyrophoric_or_highly_flammable_review", "known process gas rule", self.name, confidence="C"))

        if "C" in elems and "N" in elems:
            out.append(ClassificationRecord(chemical.chemical_id, "reactivity_flag", "combustible_or_flammable_review", "organic nitrogen screening rule", self.name, confidence="D"))
        return out

    def _matches(self, row: dict, chemical: ChemicalIdentity) -> bool:
        if row.get("inchikey") and chemical.inchikey and row["inchikey"] == chemical.inchikey:
            return True
        if row.get("cas") and chemical.cas and row["cas"] == chemical.cas:
            return True
        if row.get("name") and chemical.preferred_name and row["name"].lower() == chemical.preferred_name.lower():
            return True
        return False


def _split(text: str | None) -> list[str]:
    if not text:
        return []
    return [x.strip() for x in text.replace(";", ",").split(",") if x.strip()]
