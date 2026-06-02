from __future__ import annotations

import csv
from pathlib import Path

from gas_screening_mvp.domain.models import ChemicalIdentity, ClassificationRecord
from gas_screening_mvp.classification.pfas_rules import (
    can_evaluate_structure_rules,
    detect_fully_fluorinated_cf2_or_cf3,
    has_c_f_bond,
)


class PfasLocalClassifier:
    name = "PfasLocalClassifier"

    def __init__(self, list_csv: str | Path | None = None):
        self.list_hits: dict[str, list[str]] = {}
        if list_csv and Path(list_csv).exists():
            with Path(list_csv).open(newline="", encoding="utf-8-sig") as f:
                for row in csv.DictReader(f):
                    ik = row.get("inchikey")
                    source_list = row.get("list_name") or row.get("source") or "local_pfas_list"
                    if ik:
                        self.list_hits.setdefault(ik, []).append(source_list)

    def classify(self, chemical: ChemicalIdentity) -> list[ClassificationRecord]:
        smiles = chemical.canonical_smiles or chemical.isomeric_smiles
        rules_available = can_evaluate_structure_rules(smiles)
        motif_hit, motifs = detect_fully_fluorinated_cf2_or_cf3(smiles)
        cf = has_c_f_bond(smiles)
        hits = sorted(set(self.list_hits.get(chemical.inchikey or "", [])))

        records: list[ClassificationRecord] = []
        if motif_hit or hits:
            flag = "yes"
            basis = "structure_rule_and_list_hit" if motif_hit and hits else ("structure_rule" if motif_hit else "list_hit")
            conf = "A" if motif_hit and hits else "B"
        elif rules_available and cf:
            flag = "possible"
            basis = "contains_C_F_bond_but_no_CF2_CF3_screening_motif"
            conf = "C"
        elif rules_available and smiles:
            flag = "no"
            basis = "structure_rule_negative_and_no_local_list_hit"
            conf = "B"
        elif smiles:
            flag = "unknown"
            basis = "structure_rule_unavailable_or_unparseable"
            conf = "D"
        else:
            flag = "unknown"
            basis = "structure_unresolved"
            conf = "D"

        records.append(ClassificationRecord(chemical.chemical_id, "pfas_flag", flag, basis, self.name, confidence=conf))
        records.append(ClassificationRecord(chemical.chemical_id, "pfas_basis", basis, basis, self.name, confidence=conf))
        records.append(ClassificationRecord(chemical.chemical_id, "pfas_definition", "OECD_2021_screening_CF2_CF3_and_SEMI_CF2_CF3", basis, self.name, confidence="C"))
        for motif in motifs:
            records.append(ClassificationRecord(chemical.chemical_id, "pfas_motif", motif, "RDKit structural motif", self.name, confidence="B"))
        for hit in hits:
            records.append(ClassificationRecord(chemical.chemical_id, "pfas_list_hit", hit, "local list InChIKey match", self.name, confidence="A"))
        return records
