from __future__ import annotations

from dataclasses import dataclass, field
from collections import Counter
import re


ATOMIC_WEIGHTS = {
    "H": 1.008,
    "He": 4.0026,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.180,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "Ar": 39.948,
    "Ti": 47.867,
    "Br": 79.904,
    "As": 74.922,
    "Mo": 95.95,
    "W": 183.84,
}

DEFAULT_VALENCES = {
    "B": 3,
    "C": 4,
    "N": 3,
    "O": 2,
    "P": 3,
    "S": 2,
    "Si": 4,
    "As": 3,
}

TWO_LETTER_SYMBOLS = {
    "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Cl", "Ar",
    "Ca", "Ti", "Fe", "Br", "As", "Mo",
}

AROMATIC_SYMBOLS = {
    "b": "B",
    "c": "C",
    "n": "N",
    "o": "O",
    "p": "P",
    "s": "S",
}


@dataclass(frozen=True)
class LightweightAtom:
    index: int
    symbol: str
    explicit_hydrogens: int = 0
    bracketed: bool = False
    aromatic: bool = False


@dataclass
class LightweightSmiles:
    atoms: list[LightweightAtom] = field(default_factory=list)
    bonds: list[tuple[int, int, float]] = field(default_factory=list)
    unsupported_tokens: list[str] = field(default_factory=list)

    @property
    def element_symbols(self) -> tuple[str, ...]:
        return tuple(sorted({atom.symbol for atom in self.atoms}))

    @property
    def heavy_atom_count(self) -> int:
        return sum(1 for atom in self.atoms if atom.symbol != "H")

    def neighbors(self, atom_index: int) -> list[int]:
        out: list[int] = []
        for a, b, _order in self.bonds:
            if a == atom_index:
                out.append(b)
            elif b == atom_index:
                out.append(a)
        return out

    def bond_order_sum(self, atom_index: int) -> float:
        total = 0.0
        for a, b, order in self.bonds:
            if a == atom_index or b == atom_index:
                total += order
        return total

    def implicit_hydrogens(self, atom: LightweightAtom) -> int:
        if atom.bracketed or atom.symbol == "H":
            return 0
        valence = DEFAULT_VALENCES.get(atom.symbol)
        if valence is None:
            return 0
        if atom.aromatic and atom.symbol == "C":
            valence = 3
        used = self.bond_order_sum(atom.index)
        return max(0, int(round(valence - used)))

    def total_hydrogens(self, atom: LightweightAtom) -> int:
        return atom.explicit_hydrogens + self.implicit_hydrogens(atom)

    @property
    def formula_counts(self) -> Counter[str]:
        counts: Counter[str] = Counter()
        for atom in self.atoms:
            counts[atom.symbol] += 1
            h_count = self.total_hydrogens(atom)
            if h_count:
                counts["H"] += h_count
        return counts

    @property
    def formula(self) -> str:
        return format_formula(self.formula_counts)

    @property
    def molecular_weight(self) -> float | None:
        total = 0.0
        for symbol, count in self.formula_counts.items():
            weight = ATOMIC_WEIGHTS.get(symbol)
            if weight is None:
                return None
            total += weight * count
        return total


def parse_smiles_lightweight(smiles: str | None) -> LightweightSmiles | None:
    """Parse enough SMILES for offline MVP screening when RDKit is unavailable.

    The parser intentionally supports only atom identity, simple branching,
    ring closures, and bond order. It is not a replacement for RDKit
    canonicalization or stereochemistry handling.
    """
    if not smiles:
        return None
    text = smiles.strip()
    if not text:
        return None

    parsed = LightweightSmiles()
    branch_stack: list[int | None] = []
    ring_open: dict[str, tuple[int, float]] = {}
    current_atom: int | None = None
    pending_bond_order = 1.0
    i = 0

    while i < len(text):
        ch = text[i]
        if ch in "-:/\\":
            pending_bond_order = 1.0
            i += 1
            continue
        if ch == "=":
            pending_bond_order = 2.0
            i += 1
            continue
        if ch == "#":
            pending_bond_order = 3.0
            i += 1
            continue
        if ch == ".":
            current_atom = None
            pending_bond_order = 1.0
            i += 1
            continue
        if ch == "(":
            branch_stack.append(current_atom)
            i += 1
            continue
        if ch == ")":
            current_atom = branch_stack.pop() if branch_stack else None
            i += 1
            continue
        if ch.isdigit():
            if current_atom is None:
                parsed.unsupported_tokens.append(ch)
            elif ch in ring_open:
                other_atom, order = ring_open.pop(ch)
                parsed.bonds.append((other_atom, current_atom, pending_bond_order or order))
            else:
                ring_open[ch] = (current_atom, pending_bond_order)
            pending_bond_order = 1.0
            i += 1
            continue
        if ch == "[":
            end = text.find("]", i + 1)
            if end == -1:
                return None
            atom = _parse_bracket_atom(text[i + 1:end], len(parsed.atoms))
            if atom is None:
                parsed.unsupported_tokens.append(text[i:end + 1])
            else:
                parsed.atoms.append(atom)
                if current_atom is not None:
                    parsed.bonds.append((current_atom, atom.index, pending_bond_order))
                current_atom = atom.index
            pending_bond_order = 1.0
            i = end + 1
            continue

        atom, consumed = _parse_simple_atom(text, i, len(parsed.atoms))
        if atom is None:
            parsed.unsupported_tokens.append(ch)
            i += 1
            continue
        parsed.atoms.append(atom)
        if current_atom is not None:
            parsed.bonds.append((current_atom, atom.index, pending_bond_order))
        current_atom = atom.index
        pending_bond_order = 1.0
        i += consumed

    if not parsed.atoms or ring_open:
        return None
    return parsed


def format_formula(counts: Counter[str]) -> str:
    if not counts:
        return ""
    order: list[str] = []
    if "C" in counts:
        order.append("C")
        if "H" in counts:
            order.append("H")
    for symbol in sorted(counts):
        if symbol not in order:
            order.append(symbol)
    return "".join(symbol + (str(counts[symbol]) if counts[symbol] != 1 else "") for symbol in order)


def _parse_simple_atom(text: str, pos: int, index: int) -> tuple[LightweightAtom | None, int]:
    two = text[pos:pos + 2]
    if two in TWO_LETTER_SYMBOLS:
        return LightweightAtom(index=index, symbol=two), 2
    ch = text[pos]
    if ch in AROMATIC_SYMBOLS:
        return LightweightAtom(index=index, symbol=AROMATIC_SYMBOLS[ch], aromatic=True), 1
    if ch.isalpha() and ch.upper() == ch:
        return LightweightAtom(index=index, symbol=ch), 1
    return None, 0


def _parse_bracket_atom(content: str, index: int) -> LightweightAtom | None:
    body = re.sub(r"^\d+", "", content.strip())
    if not body:
        return None
    match = re.match(r"([A-Z][a-z]?|[bcnops])", body)
    if not match:
        return None
    raw_symbol = match.group(1)
    symbol = AROMATIC_SYMBOLS.get(raw_symbol, raw_symbol)
    h_match = re.search(r"H(\d*)", body[match.end():])
    explicit_h = 0
    if h_match:
        explicit_h = int(h_match.group(1) or "1")
    return LightweightAtom(
        index=index,
        symbol=symbol,
        explicit_hydrogens=explicit_h,
        bracketed=True,
        aromatic=raw_symbol in AROMATIC_SYMBOLS,
    )
