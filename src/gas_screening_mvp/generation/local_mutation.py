from __future__ import annotations

from typing import Iterable
import uuid

from gas_screening_mvp.domain.models import MoleculeCandidate
from gas_screening_mvp.normalization.lightweight_smiles import parse_smiles_lightweight

try:
    from rdkit import Chem
    from rdkit.Chem import inchi as rd_inchi
except Exception:  # pragma: no cover
    Chem = None
    rd_inchi = None


class LocalMutationGenerator:
    """One-generation local mutations from seed molecules.

    RDKit is used for sanitize and InChIKey dedupe when available. The fallback
    parser keeps offline tests deterministic but does not replace RDKit for
    production exploration runs.
    """

    name = "LocalMutationGenerator"

    def __init__(
        self,
        seeds: list[MoleculeCandidate],
        enabled: bool = False,
        max_candidates: int = 100,
        max_heavy_atoms: int = 20,
        allowed_families: list[str] | None = None,
    ):
        self.seeds = seeds
        self.enabled = enabled
        self.max_candidates = int(max_candidates)
        self.max_heavy_atoms = int(max_heavy_atoms)
        self.allowed_families = set(allowed_families or ["amine", "fluorocarbon"])

    def generate(self) -> Iterable[MoleculeCandidate]:
        if not self.enabled:
            return []

        out: list[MoleculeCandidate] = []
        seen: set[str] = set()
        for seed in self.seeds:
            if len(out) >= self.max_candidates:
                break
            if seed.family not in self.allowed_families or not seed.smiles:
                continue
            for rule, smiles in self._mutations(seed):
                if len(out) >= self.max_candidates:
                    break
                normalized = _sanitize_and_key(smiles)
                if normalized is None:
                    continue
                clean_smiles, dedupe_key = normalized
                if dedupe_key in seen or not self._within_limits(clean_smiles):
                    continue
                seen.add(dedupe_key)
                out.append(self._candidate(seed, rule, clean_smiles))
        return out

    def _mutations(self, seed: MoleculeCandidate) -> list[tuple[str, str]]:
        smiles = seed.smiles or ""
        if seed.family == "amine":
            return _amine_mutations(smiles)
        if seed.family == "fluorocarbon":
            return _fluorocarbon_mutations(smiles)
        return []

    def _within_limits(self, smiles: str) -> bool:
        parsed = parse_smiles_lightweight(smiles)
        return parsed is not None and not parsed.unsupported_tokens and parsed.heavy_atom_count <= self.max_heavy_atoms

    def _candidate(self, seed: MoleculeCandidate, rule: str, smiles: str) -> MoleculeCandidate:
        cid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"mutation:{seed.candidate_id}:{rule}:{smiles}"))
        return MoleculeCandidate(
            candidate_id=cid,
            source="local_mutation",
            family=seed.family,
            input_name=f"{seed.input_name or seed.candidate_id}_{rule}",
            smiles=smiles,
            cas=None,
            generation_rule=rule,
            parent_candidate_id=seed.candidate_id,
            generation_score=max(0.0, seed.generation_score - 0.15),
        )


def _amine_mutations(smiles: str) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []
    replacements = [
        ("CN", "FCN", "H_to_F_methylamine"),
        ("CN", "C(F)(F)(F)N", "CH3_to_CF3"),
        ("CCN", "FCCN", "H_to_F_ethylamine"),
        ("CCN", "FC(F)(F)CN", "ethyl_to_trifluoroethyl"),
        ("N(C)C", "N(C)CF", "amine_low_molecular_substitution"),
        ("N(C)C", "N(C)C(F)(F)F", "amine_low_molecular_substitution_cf3"),
        ("CCN(CC)CC", "CCN(CC)C", "amine_family_small_substitution"),
    ]
    for old, new, rule in replacements:
        if old in smiles:
            candidates.append((rule, smiles.replace(old, new, 1)))
    if smiles in {"CN", "CCN", "CNC", "CN(C)C"}:
        candidates.extend([
            ("amine_family_small_substitution", "CCN"),
            ("amine_family_small_substitution", "CNC"),
            ("amine_family_small_substitution", "CN(C)C"),
        ])
    return candidates


def _fluorocarbon_mutations(smiles: str) -> list[tuple[str, str]]:
    replacements = [
        ("FC(F)(F)F", "FC(F)F", "fluorocarbon_small_substitution"),
        ("FC(F)F", "FCF", "fluorocarbon_small_substitution"),
        ("FCF", "CF", "fluorocarbon_small_substitution"),
        ("CF", "FC(F)F", "H_to_F"),
        ("FC(F)(F)C(F)(F)F", "FC(F)(F)C(F)F", "fluorocarbon_small_substitution"),
        ("FC(F)(F)C(F)F", "FC(F)C(F)F", "fluorocarbon_small_substitution"),
    ]
    return [(rule, smiles.replace(old, new, 1)) for old, new, rule in replacements if old in smiles]


def _sanitize_and_key(smiles: str) -> tuple[str, str] | None:
    if Chem is not None:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return None
        clean = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        key = clean
        if rd_inchi is not None:
            try:
                inchi = rd_inchi.MolToInchi(mol)
                key = rd_inchi.InchiToInchiKey(inchi) if inchi else clean
            except Exception:
                key = clean
        return clean, key

    parsed = parse_smiles_lightweight(smiles)
    if parsed is None or parsed.unsupported_tokens:
        return None
    return smiles, smiles
