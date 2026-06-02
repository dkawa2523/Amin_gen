from __future__ import annotations

from pathlib import Path
from typing import Iterable
import csv
import uuid

from gas_screening_mvp.domain.models import MoleculeCandidate


DEFAULT_SEEDS = [
    # Amines / organic N compounds
    {"name": "methylamine", "smiles": "CN", "family": "amine", "cas": "74-89-5"},
    {"name": "dimethylamine", "smiles": "CNC", "family": "amine", "cas": "124-40-3"},
    {"name": "trimethylamine", "smiles": "CN(C)C", "family": "amine", "cas": "75-50-3"},
    {"name": "ethylamine", "smiles": "CCN", "family": "amine", "cas": "75-04-7"},
    {"name": "diethylamine", "smiles": "CCNCC", "family": "amine", "cas": "109-89-7"},
    {"name": "triethylamine", "smiles": "CCN(CC)CC", "family": "amine", "cas": "121-44-8"},
    {"name": "tert-butylamine", "smiles": "CC(C)(C)N", "family": "amine", "cas": "75-64-9"},
    {"name": "pyridine", "smiles": "c1ccncc1", "family": "amine", "cas": "110-86-1"},
    {"name": "morpholine", "smiles": "C1COCCN1", "family": "amine", "cas": "110-91-8"},
    {"name": "piperidine", "smiles": "C1CCNCC1", "family": "amine", "cas": "110-89-4"},
    # Inorganic / process gases
    {"name": "ammonia", "smiles": "N", "family": "inorganic", "cas": "7664-41-7"},
    {"name": "hydrogen fluoride", "smiles": "[H]F", "family": "inorganic", "cas": "7664-39-3"},
    {"name": "hydrogen chloride", "smiles": "[H]Cl", "family": "inorganic", "cas": "7647-01-0"},
    {"name": "chlorine", "smiles": "ClCl", "family": "inorganic", "cas": "7782-50-5"},
    {"name": "fluorine", "smiles": "FF", "family": "inorganic", "cas": "7782-41-4"},
    {"name": "boron trichloride", "smiles": "ClB(Cl)Cl", "family": "inorganic", "cas": "10294-34-5"},
    {"name": "boron trifluoride", "smiles": "FB(F)F", "family": "inorganic", "cas": "7637-07-2"},
    {"name": "silane", "smiles": "[SiH4]", "family": "inorganic", "cas": "7803-62-5"},
    {"name": "disilane", "smiles": "[SiH3][SiH3]", "family": "inorganic", "cas": "1590-87-0"},
    {"name": "phosphine", "smiles": "P", "family": "inorganic", "cas": "7803-51-2"},
    {"name": "arsine", "smiles": "[AsH3]", "family": "inorganic", "cas": "7784-42-1"},
    {"name": "tungsten hexafluoride", "smiles": "F[W](F)(F)(F)(F)F", "family": "inorganic", "cas": "7783-82-6"},
    {"name": "nitrogen trifluoride", "smiles": "FN(F)F", "family": "inorganic", "cas": "7783-54-2"},
    {"name": "nitrous oxide", "smiles": "[N-]=[N+]=O", "family": "inorganic", "cas": "10024-97-2"},
    {"name": "carbon monoxide", "smiles": "[C-]#[O+]", "family": "inorganic", "cas": "630-08-0"},
    {"name": "carbon dioxide", "smiles": "O=C=O", "family": "inorganic", "cas": "124-38-9"},
    {"name": "nitrogen", "smiles": "N#N", "family": "inorganic", "cas": "7727-37-9"},
    {"name": "oxygen", "smiles": "O=O", "family": "inorganic", "cas": "7782-44-7"},
    {"name": "argon", "smiles": "[Ar]", "family": "inorganic", "cas": "7440-37-1"},
    {"name": "helium", "smiles": "[He]", "family": "inorganic", "cas": "7440-59-7"},
    {"name": "hydrogen", "smiles": "[H][H]", "family": "inorganic", "cas": "1333-74-0"},
    # Fluorocarbons / etch-cleaning gases
    {"name": "carbon tetrafluoride", "smiles": "FC(F)(F)F", "family": "fluorocarbon", "cas": "75-73-0"},
    {"name": "trifluoromethane", "smiles": "FC(F)F", "family": "fluorocarbon", "cas": "75-46-7"},
    {"name": "difluoromethane", "smiles": "FCF", "family": "fluorocarbon", "cas": "75-10-5"},
    {"name": "fluoromethane", "smiles": "CF", "family": "fluorocarbon", "cas": "593-53-3"},
    {"name": "hexafluoroethane", "smiles": "FC(F)(F)C(F)(F)F", "family": "fluorocarbon", "cas": "76-16-4"},
    {"name": "octafluoropropane", "smiles": "FC(F)(F)C(F)(F)C(F)(F)F", "family": "fluorocarbon", "cas": "76-19-7"},
    {"name": "octafluorocyclobutane", "smiles": "C1(F)(F)C(F)(F)C(F)(F)C1(F)F", "family": "fluorocarbon", "cas": "115-25-3"},
    {"name": "sulfur hexafluoride", "smiles": "FS(F)(F)(F)(F)F", "family": "fluorocarbon", "cas": "2551-62-4"},
]


class SeedListGenerator:
    name = "SeedListGenerator"

    def __init__(self, csv_path: str | Path | None = None):
        self.csv_path = Path(csv_path) if csv_path else None

    def generate(self) -> Iterable[MoleculeCandidate]:
        rows = DEFAULT_SEEDS
        if self.csv_path and self.csv_path.exists():
            with self.csv_path.open(newline="", encoding="utf-8-sig") as f:
                rows = list(csv.DictReader(f))

        for row in rows:
            name = row.get("name") or row.get("input_name") or row.get("preferred_name")
            smiles = row.get("smiles") or row.get("canonical_smiles")
            formula = row.get("formula") or None
            family = row.get("family") or "manual"
            cas = row.get("cas") or None
            cid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"seed:{name}:{smiles}:{formula}:{cas}"))
            yield MoleculeCandidate(
                candidate_id=cid,
                source="seed",
                family=family,
                input_name=name,
                smiles=smiles,
                formula=formula,
                cas=cas,
                generation_rule="seed_list",
                generation_score=1.0,
            )
