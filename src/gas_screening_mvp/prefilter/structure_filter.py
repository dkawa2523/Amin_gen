from __future__ import annotations

from dataclasses import dataclass
from gas_screening_mvp.domain.models import NormalizedMolecule
from gas_screening_mvp.prefilter.score import api_priority_score


@dataclass(frozen=True)
class FilterResult:
    molecule: NormalizedMolecule
    passed: bool
    reasons: tuple[str, ...]
    api_priority_score: float


class CandidatePrefilter:
    def __init__(self, config: dict):
        self.max_molecular_weight = float(config.get("max_molecular_weight", 300))
        self.max_heavy_atoms = int(config.get("max_heavy_atoms", 20))
        self.allow_formula_only = bool(config.get("allow_formula_only", False))
        self.min_score = float(config.get("min_score", 0.0))

    def evaluate(self, mol: NormalizedMolecule) -> FilterResult:
        reasons: list[str] = []

        if mol.structure_status != "valid":
            if not (self.allow_formula_only and mol.structure_status == "manual_review_required"):
                reasons.append(f"structure_status={mol.structure_status}")

        if mol.molecular_weight is not None and mol.molecular_weight > self.max_molecular_weight:
            reasons.append("molecular_weight_too_high")

        if mol.heavy_atom_count is not None and mol.heavy_atom_count > self.max_heavy_atoms:
            reasons.append("too_many_heavy_atoms")

        score = api_priority_score(mol)
        if score < self.min_score:
            reasons.append("priority_score_too_low")

        return FilterResult(
            molecule=mol,
            passed=not reasons,
            reasons=tuple(reasons),
            api_priority_score=score,
        )

    def filter(self, molecules: list[NormalizedMolecule]) -> tuple[list[FilterResult], list[FilterResult]]:
        passed: list[FilterResult] = []
        rejected: list[FilterResult] = []
        for mol in molecules:
            result = self.evaluate(mol)
            (passed if result.passed else rejected).append(result)
        return passed, rejected
