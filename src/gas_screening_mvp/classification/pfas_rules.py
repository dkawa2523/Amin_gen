from __future__ import annotations

try:
    from rdkit import Chem
except Exception:  # pragma: no cover
    Chem = None

from gas_screening_mvp.normalization.lightweight_smiles import parse_smiles_lightweight


def can_evaluate_structure_rules(smiles: str | None) -> bool:
    if not smiles:
        return False
    if Chem is not None:
        return Chem.MolFromSmiles(smiles) is not None
    parsed = parse_smiles_lightweight(smiles)
    return parsed is not None and not parsed.unsupported_tokens


def has_c_f_bond(smiles: str | None) -> bool:
    if not smiles:
        return False
    if Chem is not None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        for bond in mol.GetBonds():
            atoms = {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
            if atoms == {6, 9}:
                return True
        return False

    parsed = parse_smiles_lightweight(smiles)
    if parsed is None or parsed.unsupported_tokens:
        return False
    for a, b, _order in parsed.bonds:
        atoms = {parsed.atoms[a].symbol, parsed.atoms[b].symbol}
        if atoms == {"C", "F"}:
            return True
    return False


def detect_fully_fluorinated_cf2_or_cf3(smiles: str | None) -> tuple[bool, list[str]]:
    """Screen for OECD/SEMI-like -CF2- / -CF3 motifs.

    This is a screening rule, not a legal/regulatory conclusion. It flags carbon
    atoms with no attached hydrogens and at least two directly attached fluorines.
    """
    if not smiles:
        return False, []
    motifs: list[str] = []
    if Chem is not None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            total_h = atom.GetTotalNumHs()
            f_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 9)
            if total_h == 0 and f_neighbors >= 3:
                motifs.append("CF3_or_perfluorinated_C")
            elif total_h == 0 and f_neighbors == 2:
                motifs.append("CF2")
        return bool(motifs), sorted(set(motifs))

    parsed = parse_smiles_lightweight(smiles)
    if parsed is None or parsed.unsupported_tokens:
        return False, []
    for atom in parsed.atoms:
        if atom.symbol != "C":
            continue
        total_h = parsed.total_hydrogens(atom)
        f_neighbors = sum(1 for n in parsed.neighbors(atom.index) if parsed.atoms[n].symbol == "F")
        if total_h == 0 and f_neighbors >= 3:
            motifs.append("CF3_or_perfluorinated_C")
        elif total_h == 0 and f_neighbors == 2:
            motifs.append("CF2")
    return bool(motifs), sorted(set(motifs))
