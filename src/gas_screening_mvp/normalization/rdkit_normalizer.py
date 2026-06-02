from __future__ import annotations

from typing import Iterable

from gas_screening_mvp.domain.models import MoleculeCandidate, NormalizedMolecule
from gas_screening_mvp.normalization.lightweight_smiles import parse_smiles_lightweight

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.Chem import inchi as rd_inchi
except Exception:  # pragma: no cover - allows import without RDKit
    Chem = None
    Descriptors = None
    rdMolDescriptors = None
    rd_inchi = None


_ALLOWED_DEFAULT = {
    "H", "C", "N", "O", "F", "Cl", "Br", "Si", "B", "P", "S",
    "W", "Mo", "Ti", "Ar", "He", "Ne", "As"
}


class RDKitNormalizer:
    name = "RDKitNormalizer"

    def __init__(self, allowed_elements: set[str] | None = None):
        self.allowed_elements = allowed_elements or set(_ALLOWED_DEFAULT)

    def _normalize_without_rdkit(self, candidate: MoleculeCandidate) -> NormalizedMolecule:
        if not candidate.smiles:
            return _formula_only(candidate)

        smiles = candidate.smiles.strip()
        parsed = parse_smiles_lightweight(smiles)
        if parsed is None or parsed.unsupported_tokens:
            reason = "RDKit is not installed and lightweight SMILES parser could not safely parse structure"
            if parsed and parsed.unsupported_tokens:
                reason = f"{reason}: {','.join(sorted(set(parsed.unsupported_tokens)))}"
            return _invalid(candidate, reason)

        elements = parsed.element_symbols
        unsupported = [e for e in elements if e not in self.allowed_elements]
        if unsupported:
            return _unsupported(
                candidate,
                smiles,
                elements,
                f"unsupported elements: {','.join(unsupported)}",
            )

        molecular_weight = parsed.molecular_weight
        return NormalizedMolecule(
            candidate_id=candidate.candidate_id,
            source=candidate.source,
            family=candidate.family,
            input_name=candidate.input_name,
            cas=candidate.cas,
            canonical_smiles=smiles,
            isomeric_smiles=smiles,
            standard_inchi=None,
            standard_inchikey=None,
            formula=candidate.formula or parsed.formula,
            molecular_weight=float(molecular_weight) if molecular_weight is not None else None,
            heavy_atom_count=parsed.heavy_atom_count,
            element_symbols=elements,
            structure_status="valid",
            status_reason="RDKit is not installed; normalized with lightweight SMILES parser",
            metadata=_candidate_metadata(candidate),
        )

    def normalize(self, candidate: MoleculeCandidate) -> NormalizedMolecule:
        if Chem is None:
            return self._normalize_without_rdkit(candidate)

        if not candidate.smiles:
            return NormalizedMolecule(
                candidate_id=candidate.candidate_id,
                source=candidate.source,
                family=candidate.family,
                input_name=candidate.input_name,
                cas=candidate.cas,
                canonical_smiles=None,
                isomeric_smiles=None,
                standard_inchi=None,
                standard_inchikey=None,
                formula=candidate.formula,
                molecular_weight=None,
                heavy_atom_count=None,
                element_symbols=(),
                structure_status="manual_review_required",
                status_reason="formula-only candidate; structure lookup required",
                metadata=_candidate_metadata(candidate),
            )

        smiles = candidate.smiles.strip()
        if "." in smiles:
            status = "mixture"
        else:
            status = "valid"

        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return NormalizedMolecule(
                candidate_id=candidate.candidate_id,
                source=candidate.source,
                family=candidate.family,
                input_name=candidate.input_name,
                cas=candidate.cas,
                canonical_smiles=None,
                isomeric_smiles=None,
                standard_inchi=None,
                standard_inchikey=None,
                formula=candidate.formula,
                molecular_weight=None,
                heavy_atom_count=None,
                element_symbols=(),
                structure_status="invalid",
                status_reason=f"RDKit failed to parse SMILES: {smiles}",
                metadata=_candidate_metadata(candidate),
            )

        radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        if radicals:
            status = "radical"

        elements = tuple(sorted({atom.GetSymbol() for atom in mol.GetAtoms()}))
        unsupported = [e for e in elements if e not in self.allowed_elements]
        if unsupported:
            status = "unsupported"
            reason = f"unsupported elements: {','.join(unsupported)}"
        else:
            reason = None

        try:
            can = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
            iso = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        except Exception:
            can = smiles
            iso = smiles

        formula = None
        mw = None
        heavy = None
        try:
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mw = float(Descriptors.MolWt(mol))
            heavy = int(mol.GetNumHeavyAtoms())
        except Exception:
            formula = candidate.formula

        inchi = None
        ikey = None
        try:
            inchi = rd_inchi.MolToInchi(mol)
            ikey = rd_inchi.InchiToInchiKey(inchi) if inchi else None
        except Exception:
            pass

        return NormalizedMolecule(
            candidate_id=candidate.candidate_id,
            source=candidate.source,
            family=candidate.family,
            input_name=candidate.input_name,
            cas=candidate.cas,
            canonical_smiles=can,
            isomeric_smiles=iso,
            standard_inchi=inchi,
            standard_inchikey=ikey,
            formula=formula or candidate.formula,
            molecular_weight=mw,
            heavy_atom_count=heavy,
            element_symbols=elements,
            structure_status=status,
            status_reason=reason,
            metadata=_candidate_metadata(candidate),
        )


def normalize_all(candidates: Iterable[MoleculeCandidate], allowed_elements: set[str] | None = None) -> list[NormalizedMolecule]:
    norm = RDKitNormalizer(allowed_elements=allowed_elements)
    return [norm.normalize(c) for c in candidates]


def _formula_only(candidate: MoleculeCandidate) -> NormalizedMolecule:
    return NormalizedMolecule(
        candidate_id=candidate.candidate_id,
        source=candidate.source,
        family=candidate.family,
        input_name=candidate.input_name,
        cas=candidate.cas,
        canonical_smiles=None,
        isomeric_smiles=None,
        standard_inchi=None,
        standard_inchikey=None,
        formula=candidate.formula,
        molecular_weight=None,
        heavy_atom_count=None,
        element_symbols=(),
        structure_status="manual_review_required",
        status_reason="formula-only candidate; structure lookup required",
        metadata=_candidate_metadata(candidate),
    )


def _invalid(candidate: MoleculeCandidate, reason: str) -> NormalizedMolecule:
    return NormalizedMolecule(
        candidate_id=candidate.candidate_id,
        source=candidate.source,
        family=candidate.family,
        input_name=candidate.input_name,
        cas=candidate.cas,
        canonical_smiles=None,
        isomeric_smiles=None,
        standard_inchi=None,
        standard_inchikey=None,
        formula=candidate.formula,
        molecular_weight=None,
        heavy_atom_count=None,
        element_symbols=(),
        structure_status="invalid",
        status_reason=reason,
        metadata=_candidate_metadata(candidate),
    )


def _unsupported(candidate: MoleculeCandidate, smiles: str, elements: tuple[str, ...], reason: str) -> NormalizedMolecule:
    return NormalizedMolecule(
        candidate_id=candidate.candidate_id,
        source=candidate.source,
        family=candidate.family,
        input_name=candidate.input_name,
        cas=candidate.cas,
        canonical_smiles=smiles,
        isomeric_smiles=smiles,
        standard_inchi=None,
        standard_inchikey=None,
        formula=candidate.formula,
        molecular_weight=None,
        heavy_atom_count=None,
        element_symbols=elements,
        structure_status="unsupported",
        status_reason=reason,
        metadata=_candidate_metadata(candidate),
    )


def _candidate_metadata(candidate: MoleculeCandidate) -> dict[str, str]:
    metadata = dict(candidate.metadata)
    metadata.setdefault("generation_rule", candidate.generation_rule)
    if candidate.parent_candidate_id:
        metadata.setdefault("parent_candidate_id", candidate.parent_candidate_id)
    return metadata
