from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement
from typing import Iterable
import uuid

from gas_screening_mvp.domain.models import MoleculeCandidate
from gas_screening_mvp.normalization.lightweight_smiles import parse_smiles_lightweight


@dataclass(frozen=True)
class Substituent:
    label: str
    smiles_prefix: str
    carbons: int
    category: str


AMINE_CLASS_BY_RULE = {
    "primary_amine": "primary",
    "secondary_amine": "secondary",
    "tertiary_amine": "tertiary",
    "cyclic_amine": "cyclic",
}

DEFAULT_SUBSTITUENTS = [
    Substituent("methyl", "C", 1, "small_alkyl"),
    Substituent("ethyl", "CC", 2, "small_alkyl"),
    Substituent("n-propyl", "CCC", 3, "small_alkyl"),
    Substituent("isopropyl", "C(C)C", 3, "small_alkyl"),
    Substituent("n-butyl", "CCCC", 4, "alkyl"),
    Substituent("tert-butyl", "C(C)(C)C", 4, "alkyl"),
    Substituent("vinyl", "C=C", 2, "unsaturated_alkyl"),
    Substituent("allyl", "CC=C", 3, "unsaturated_alkyl"),
    Substituent("cyclopropyl", "C1CC1", 3, "cyclic_alkyl"),
    Substituent("fluoromethyl", "FC", 1, "fluorinated_alkyl"),
    Substituent("difluoromethyl", "FC(F)", 1, "fluorinated_alkyl"),
    Substituent("fluoroethyl", "FCC", 2, "fluorinated_alkyl"),
    Substituent("difluoroethyl", "FC(F)C", 2, "fluorinated_alkyl"),
    Substituent("trifluoromethyl", "FC(F)(F)", 1, "fluorinated_alkyl"),
    Substituent("trifluoroethyl", "FC(F)(F)C", 2, "fluorinated_alkyl"),
    Substituent("pentafluoroethyl", "FC(F)(F)C(F)(F)", 2, "fluorinated_alkyl"),
]

CYCLIC_AMINES = {
    "aziridine": "C1CN1",
    "azetidine": "C1CCN1",
    "pyrrolidine": "C1CCNC1",
    "piperidine": "C1CCNCC1",
    "morpholine": "C1COCCN1",
    "piperazine": "C1CNCCN1",
}


class AmineTemplateGenerator:
    name = "AmineTemplateGenerator"

    def __init__(
        self,
        max_candidates: int = 500,
        max_carbons: int = 6,
        max_heavy_atoms: int = 12,
        allowed_substituents: list[str] | None = None,
    ):
        self.max_candidates = int(max_candidates)
        self.max_carbons = int(max_carbons)
        self.max_heavy_atoms = int(max_heavy_atoms)
        allowed = {s.strip() for s in allowed_substituents or [] if s.strip()}
        self.substituents = [
            s for s in DEFAULT_SUBSTITUENTS
            if not allowed or s.label in allowed or s.category in allowed
        ]

    def generate(self) -> Iterable[MoleculeCandidate]:
        seen: set[str] = set()
        count = 0

        def emit(rule: str, smiles: str, label: str, score: float):
            nonlocal count
            if count >= self.max_candidates or smiles in seen or not self._within_limits(smiles):
                return None
            seen.add(smiles)
            count += 1
            return self._candidate(rule, smiles, label, score)

        # Primary amines: R-NH2.
        for substituent in self.substituents:
            candidate = emit("primary_amine", f"{substituent.smiles_prefix}N", substituent.label, 0.70)
            if candidate:
                yield candidate

        # Secondary amines: R1-NH-R2.
        for s1, s2 in combinations_with_replacement(self.substituents, 2):
            if s1.carbons + s2.carbons > self.max_carbons:
                continue
            candidate = emit(
                "secondary_amine",
                f"N({s1.smiles_prefix}){s2.smiles_prefix}",
                f"{s1.label}_{s2.label}",
                0.62,
            )
            if candidate:
                yield candidate

        # Tertiary amines: R1-N(R2)-R3, bounded by carbon/heavy atom filters.
        for s1, s2, s3 in combinations_with_replacement(self.substituents, 3):
            if s1.carbons + s2.carbons + s3.carbons > self.max_carbons:
                continue
            candidate = emit(
                "tertiary_amine",
                f"N({s1.smiles_prefix})({s2.smiles_prefix}){s3.smiles_prefix}",
                f"{s1.label}_{s2.label}_{s3.label}",
                0.58,
            )
            if candidate:
                yield candidate

        for label, smiles in CYCLIC_AMINES.items():
            candidate = emit("cyclic_amine", smiles, label, 0.64)
            if candidate:
                yield candidate

    def _within_limits(self, smiles: str) -> bool:
        parsed = parse_smiles_lightweight(smiles)
        if parsed is None or parsed.unsupported_tokens:
            return False
        counts = parsed.formula_counts
        if counts.get("C", 0) > self.max_carbons:
            return False
        return parsed.heavy_atom_count <= self.max_heavy_atoms

    def _candidate(self, rule: str, smiles: str, label: str, score: float) -> MoleculeCandidate:
        cid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"amine:{rule}:{smiles}"))
        amine_class = AMINE_CLASS_BY_RULE.get(rule, "unknown")
        substituents = [] if rule == "cyclic_amine" else label.split("_")
        categories = _substituent_categories(substituents)
        fluorinated_count = sum(1 for item in substituents if _is_fluorinated(item))
        metadata = {
            "amine_class": amine_class,
            "amine_substituents": "; ".join(substituents),
            "amine_substituent_count": str(len(substituents)),
            "fluorinated_substituent_count": str(fluorinated_count),
            "contains_fluorinated_substituent": str(fluorinated_count > 0).lower(),
            "contains_unsaturated_substituent": str(any(c == "unsaturated_alkyl" for c in categories)).lower(),
            "contains_cyclic_substituent": str(any(c == "cyclic_alkyl" for c in categories) or rule == "cyclic_amine").lower(),
        }
        if rule == "cyclic_amine":
            metadata["ring_name"] = label
        return MoleculeCandidate(
            candidate_id=cid,
            source="template_generation",
            family="amine",
            input_name=f"generated_{rule}_{label}",
            smiles=smiles,
            generation_rule=rule,
            generation_score=score,
            metadata=metadata,
        )


def _substituent_categories(labels: list[str]) -> list[str]:
    by_label = {s.label: s.category for s in DEFAULT_SUBSTITUENTS}
    return [by_label.get(label, "") for label in labels]


def _is_fluorinated(label: str) -> bool:
    category = {s.label: s.category for s in DEFAULT_SUBSTITUENTS}.get(label, "")
    return category == "fluorinated_alkyl" or "fluoro" in label
