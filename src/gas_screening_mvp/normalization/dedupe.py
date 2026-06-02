from __future__ import annotations

from gas_screening_mvp.domain.models import NormalizedMolecule


def dedupe_molecules(molecules: list[NormalizedMolecule]) -> list[NormalizedMolecule]:
    """Deduplicate molecules by InChIKey when available, else by structure/name.

    Keeps the highest-information record first. The generator should pass seed
    candidates before lower-confidence generated candidates when seed records are
    preferred.
    """
    seen: set[str] = set()
    out: list[NormalizedMolecule] = []
    for mol in molecules:
        key = mol.standard_inchikey or mol.canonical_smiles or f"{mol.formula}:{mol.input_name}"
        if key in seen:
            continue
        seen.add(key)
        out.append(mol)
    return out
