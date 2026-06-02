from __future__ import annotations

from typing import Iterable
import uuid

from gas_screening_mvp.domain.models import MoleculeCandidate


CURATED_FLUOROCARBON_SMILES = {
    "CF4": "FC(F)(F)F",
    "CHF3": "FC(F)F",
    "CH2F2": "FCF",
    "CH3F": "CF",
    "C2F6": "FC(F)(F)C(F)(F)F",
    "C2HF5": "FC(F)(F)C(F)F",
    "C2H2F4": "FC(F)C(F)F",
    "C3F8": "FC(F)(F)C(F)(F)C(F)(F)F",
    "C4F8_cyclic": "C1(F)(F)C(F)(F)C(F)(F)C1(F)F",
    "C4F6": "FC(F)=C(F)C(F)=C(F)F",  # screening candidate, isomer not uniquely specified
    "NF3": "FN(F)F",
    "SF6": "FS(F)(F)(F)(F)F",
}


class FluorocarbonFormulaGenerator:
    name = "FluorocarbonFormulaGenerator"

    def __init__(self, include_formula_only: bool = True, carbon_range: tuple[int, int] | list[int] | None = None):
        self.include_formula_only = include_formula_only
        if carbon_range and len(carbon_range) == 2:
            self.carbon_min = int(carbon_range[0])
            self.carbon_max = int(carbon_range[1])
        else:
            self.carbon_min = 1
            self.carbon_max = 5

    def generate(self) -> Iterable[MoleculeCandidate]:
        for label, smiles in CURATED_FLUOROCARBON_SMILES.items():
            cid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"fluoro:{label}:{smiles}"))
            yield MoleculeCandidate(
                candidate_id=cid,
                source="formula_generation",
                family="fluorocarbon",
                input_name=label,
                smiles=smiles,
                formula=label.replace("_cyclic", ""),
                generation_rule="curated_fluorocarbon_formula",
                generation_score=0.78,
                metadata={
                    "candidate_scope": "curated_fluorinated_process_gas",
                    "semiconductor_process_roles": "etch_clean_candidate; greenhouse_gas_review",
                },
            )

        if self.include_formula_only:
            # Formula-only candidates are intended for PubChem formula search in
            # exploration mode. They are not sent to structure-dependent local rules.
            for c in range(self.carbon_min, self.carbon_max + 1):
                for h in range(0, 2 * c + 2):
                    for f in range(1, 2 * c + 3):
                        formula = f"C{c if c > 1 else ''}H{h if h > 1 else ('' if h == 1 else '0')}F{f if f > 1 else ''}"
                        if "H0" in formula:
                            formula = formula.replace("H0", "")
                        cid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"fluoro_formula:{formula}"))
                        yield MoleculeCandidate(
                            candidate_id=cid,
                            source="formula_generation",
                            family="fluorocarbon",
                            input_name=formula,
                            smiles=None,
                            formula=formula,
                            generation_rule="fluorocarbon_formula_only",
                            generation_score=0.35,
                        )
