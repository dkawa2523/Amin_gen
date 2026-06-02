from __future__ import annotations

from gas_screening_mvp.domain.models import NormalizedMolecule


FAMILY_PRIORITY = {
    "inorganic": 0.25,
    "fluorocarbon": 0.22,
    "amine": 0.20,
    "organosilicon": 0.18,
    "manual": 0.15,
}

SOURCE_PRIORITY = {
    "seed": 0.35,
    "semiconductor_amine": 0.30,
    "pubchem_similarity": 0.25,
    "pubchem_formula": 0.20,
    "template_generation": 0.12,
    "formula_generation": 0.10,
    "local_mutation": 0.10,
    "manual": 0.25,
}


def api_priority_score(mol: NormalizedMolecule) -> float:
    score = 0.0
    score += FAMILY_PRIORITY.get(mol.family, 0.10)
    score += SOURCE_PRIORITY.get(mol.source, 0.05)

    if mol.molecular_weight is not None:
        if mol.molecular_weight <= 120:
            score += 0.18
        elif mol.molecular_weight <= 200:
            score += 0.10
        elif mol.molecular_weight <= 300:
            score += 0.03
        else:
            score -= 0.20

    if mol.heavy_atom_count is not None:
        if mol.heavy_atom_count <= 8:
            score += 0.12
        elif mol.heavy_atom_count <= 15:
            score += 0.06
        else:
            score -= 0.10

    elems = set(mol.element_symbols)
    if "F" in elems and "C" in elems:
        score += 0.10  # PFAS/GWP interest for fluorinated candidates
    if "F" in elems and ({"S", "N"} & elems):
        score += 0.06  # semiconductor fluorinated process-gas relevance
    if {"Si", "B", "P", "As", "W", "Mo", "Ti"} & elems:
        score += 0.08  # semiconductor process-gas relevance
    if mol.standard_inchikey:
        score += 0.05
    else:
        score -= 0.10
    if mol.source == "formula_generation" and mol.formula and not mol.canonical_smiles:
        score += 0.30
    if mol.metadata.get("pubchem_cid"):
        score += 0.55
    if mol.metadata.get("candidate_scope") == "semiconductor_amine_curated" or mol.metadata.get("semiconductor_process_roles"):
        score += 0.12
    if mol.structure_status != "valid":
        if mol.metadata.get("pubchem_cid") or (mol.source == "formula_generation" and mol.formula):
            score -= 0.10
        else:
            score -= 0.30

    return max(0.0, min(1.0, score))
