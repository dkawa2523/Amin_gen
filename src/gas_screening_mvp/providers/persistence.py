from __future__ import annotations

from gas_screening_mvp.domain.models import ChemicalIdentity, ClassificationRecord


class PersistenceScreener:
    """Screening-level persistence classifier.

    The MVP intentionally does not assign precise half-lives unless a validated
    external workflow is added. PFAS-positive structures are flagged as likely
    persistent for screening; all model-based values should remain `model_only`.
    """

    name = "PersistenceScreener"

    def __init__(self, pfas_records_by_chemical: dict[str, list[ClassificationRecord]] | None = None):
        self.pfas_records_by_chemical = pfas_records_by_chemical or {}

    def classify(self, chemical: ChemicalIdentity) -> list[ClassificationRecord]:
        pfas_flag = None
        for rec in self.pfas_records_by_chemical.get(chemical.chemical_id, []):
            if rec.classification_type == "pfas_flag":
                pfas_flag = rec.value
                break
        records: list[ClassificationRecord] = []
        if pfas_flag == "yes":
            value = "likely_persistent"
            basis = "PFAS structure/list screening positive"
            conf = "B"
        elif pfas_flag == "possible":
            value = "unknown"
            basis = "fluorinated structure requires review"
            conf = "C"
        elif chemical.formula and set(chemical.formula) <= set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789") and not chemical.canonical_smiles:
            value = "not_applicable"
            basis = "structure unavailable or inorganic/simple formula; persistence model not run"
            conf = "D"
        else:
            value = "unknown"
            basis = "no validated persistence model configured in MVP"
            conf = "D"
        records.append(ClassificationRecord(chemical.chemical_id, "persistence_screen", value, basis, self.name, confidence=conf))
        records.append(ClassificationRecord(chemical.chemical_id, "persistence_basis", basis, basis, self.name, confidence=conf))
        return records
