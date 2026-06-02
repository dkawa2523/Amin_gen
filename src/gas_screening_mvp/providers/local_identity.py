from __future__ import annotations

import uuid
from gas_screening_mvp.domain.models import NormalizedMolecule, ChemicalIdentity


class LocalIdentityProvider:
    """Creates an identity from normalized local structure.

    This provider is deliberately conservative: it does not claim that a
    generated molecule is known commercially or registered in a public DB.
    Remote PubChem/CAS providers can enrich or override CID/CAS later.
    """

    name = "LocalIdentityProvider"

    def resolve(self, molecule: NormalizedMolecule) -> list[ChemicalIdentity]:
        if molecule.structure_status not in {"valid", "manual_review_required"}:
            return [
                ChemicalIdentity(
                    chemical_id=str(uuid.uuid5(uuid.NAMESPACE_URL, f"unresolved:{molecule.candidate_id}")),
                    candidate_id=molecule.candidate_id,
                    preferred_name=molecule.input_name,
                    cas=molecule.cas,
                    pubchem_cid=None,
                    formula=molecule.formula,
                    molecular_weight=molecule.molecular_weight,
                    canonical_smiles=molecule.canonical_smiles,
                    isomeric_smiles=molecule.isomeric_smiles,
                    inchi=molecule.standard_inchi,
                    inchikey=molecule.standard_inchikey,
                    identity_status="manual_review_required",
                    confidence=0.2,
                    source=self.name,
                )
            ]

        key = molecule.standard_inchikey or molecule.canonical_smiles or molecule.candidate_id
        return [
            ChemicalIdentity(
                chemical_id=str(uuid.uuid5(uuid.NAMESPACE_URL, f"chemical:{key}")),
                candidate_id=molecule.candidate_id,
                preferred_name=molecule.input_name,
                cas=molecule.cas,
                pubchem_cid=None,
                formula=molecule.formula,
                molecular_weight=molecule.molecular_weight,
                canonical_smiles=molecule.canonical_smiles,
                isomeric_smiles=molecule.isomeric_smiles,
                inchi=molecule.standard_inchi,
                inchikey=molecule.standard_inchikey,
                identity_status="resolved" if molecule.standard_inchikey else "manual_review_required",
                confidence=0.65 if molecule.standard_inchikey else 0.35,
                source=self.name,
            )
        ]
