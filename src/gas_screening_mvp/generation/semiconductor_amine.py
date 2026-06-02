from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
import uuid

from gas_screening_mvp.domain.models import MoleculeCandidate
from gas_screening_mvp.normalization.lightweight_smiles import parse_smiles_lightweight


@dataclass(frozen=True)
class SemiconductorAmineSeed:
    name: str
    smiles: str
    amine_class: str
    substituents: tuple[str, ...]
    precursor_family: str
    process_roles: tuple[str, ...]
    relevance_basis: str
    cas: str | None = None


SEMICONDUCTOR_AMINE_SEEDS = [
    # Volatile low-molecular amines and basic nitrogen sources.
    SemiconductorAmineSeed("ammonia", "N", "inorganic_nitrogen_source", (), "nitrogen_source", ("nitridation", "ald_nitrogen_source", "chamber_clean_review"), "common semiconductor nitrogen-containing process gas", "7664-41-7"),
    SemiconductorAmineSeed("methylamine", "CN", "primary", ("methyl",), "volatile_amine", ("nitrogen_source_candidate", "volatile_base_candidate"), "small volatile primary amine", "74-89-5"),
    SemiconductorAmineSeed("dimethylamine", "CNC", "secondary", ("methyl", "methyl"), "volatile_amine", ("nitrogen_source_candidate", "volatile_base_candidate"), "small volatile secondary amine", "124-40-3"),
    SemiconductorAmineSeed("trimethylamine", "CN(C)C", "tertiary", ("methyl", "methyl", "methyl"), "volatile_amine", ("volatile_base_candidate", "surface_chemistry_modifier"), "small volatile tertiary amine", "75-50-3"),
    SemiconductorAmineSeed("ethylamine", "CCN", "primary", ("ethyl",), "volatile_amine", ("nitrogen_source_candidate", "volatile_base_candidate"), "small volatile primary amine", "75-04-7"),
    SemiconductorAmineSeed("diethylamine", "N(CC)CC", "secondary", ("ethyl", "ethyl"), "volatile_amine", ("ligand_exchange_candidate", "volatile_base_candidate"), "small volatile secondary amine", "109-89-7"),
    SemiconductorAmineSeed("triethylamine", "CCN(CC)CC", "tertiary", ("ethyl", "ethyl", "ethyl"), "volatile_amine", ("volatile_base_candidate", "surface_chemistry_modifier"), "small volatile tertiary amine", "121-44-8"),
    SemiconductorAmineSeed("isopropylamine", "CC(C)N", "primary", ("isopropyl",), "volatile_amine", ("nitrogen_source_candidate", "volatile_base_candidate"), "branched volatile primary amine", "75-31-0"),
    SemiconductorAmineSeed("tert-butylamine", "CC(C)(C)N", "primary", ("tert-butyl",), "bulky_amine", ("sterically_hindered_ligand_candidate", "volatile_base_candidate"), "bulky primary amine used as ligand motif", "75-64-9"),
    SemiconductorAmineSeed("diisopropylamine", "N(C(C)C)C(C)C", "secondary", ("isopropyl", "isopropyl"), "bulky_amine", ("sterically_hindered_ligand_candidate", "volatile_base_candidate"), "bulky secondary amine ligand motif", "108-18-9"),
    SemiconductorAmineSeed("diisopropylethylamine", "N(CC)(C(C)C)C(C)C", "tertiary", ("ethyl", "isopropyl", "isopropyl"), "bulky_amine", ("sterically_hindered_base_candidate", "ligand_exchange_candidate"), "hindered tertiary amine base motif", "7087-68-5"),
    # Diamines and chelating amines used as ligand families.
    SemiconductorAmineSeed("ethylenediamine", "NCCN", "diamine", ("aminoethyl",), "diamine_ligand", ("chelating_ligand_candidate", "nitrogen_source_candidate"), "small diamine ligand family", "107-15-3"),
    SemiconductorAmineSeed("1,3-diaminopropane", "NCCCN", "diamine", ("aminopropyl",), "diamine_ligand", ("chelating_ligand_candidate", "nitrogen_source_candidate"), "small diamine ligand family", "109-76-2"),
    SemiconductorAmineSeed("N,N-dimethylethylenediamine", "CN(C)CCN", "diamine", ("dimethylaminoethyl",), "diamine_ligand", ("chelating_ligand_candidate", "surface_chemistry_modifier"), "tertiary/primary diamine ligand motif", "108-00-9"),
    SemiconductorAmineSeed("tetramethylethylenediamine", "CN(C)CCN(C)C", "diamine", ("dimethylaminoethyl", "methyl", "methyl"), "diamine_ligand", ("chelating_ligand_candidate", "volatile_base_candidate"), "tertiary diamine ligand motif", "110-18-9"),
    SemiconductorAmineSeed("tris(2-aminoethyl)amine", "N(CCN)(CCN)CCN", "polyamine", ("aminoethyl", "aminoethyl", "aminoethyl"), "polyamine_ligand", ("multidentate_ligand_candidate", "nitrogen_source_candidate"), "multidentate amine ligand motif", "4097-89-6"),
    # Cyclic amines and saturated N-heterocycles.
    SemiconductorAmineSeed("pyrrolidine", "C1CCNC1", "cyclic", (), "cyclic_amine", ("volatile_base_candidate", "ligand_exchange_candidate"), "small saturated cyclic amine", "123-75-1"),
    SemiconductorAmineSeed("piperidine", "C1CCNCC1", "cyclic", (), "cyclic_amine", ("volatile_base_candidate", "ligand_exchange_candidate"), "small saturated cyclic amine", "110-89-4"),
    SemiconductorAmineSeed("morpholine", "C1COCCN1", "cyclic", (), "cyclic_amine", ("oxygenated_ligand_candidate", "volatile_base_candidate"), "oxygen-containing cyclic amine", "110-91-8"),
    SemiconductorAmineSeed("piperazine", "C1CNCCN1", "cyclic", (), "cyclic_amine", ("diamine_ligand_candidate", "nitrogen_source_candidate"), "cyclic diamine ligand motif", "110-85-0"),
    SemiconductorAmineSeed("N-methylpyrrolidine", "CN1CCCC1", "cyclic", ("methyl",), "cyclic_amine", ("volatile_base_candidate", "surface_chemistry_modifier"), "N-alkyl cyclic amine motif", "120-94-5"),
    SemiconductorAmineSeed("N-methylpiperidine", "CN1CCCCC1", "cyclic", ("methyl",), "cyclic_amine", ("volatile_base_candidate", "surface_chemistry_modifier"), "N-alkyl cyclic amine motif", "626-67-5"),
    # Fluorinated amines for etch/clean and PFAS review workflows.
    SemiconductorAmineSeed("2-fluoroethylamine", "FCCN", "primary", ("fluoroethyl",), "fluorinated_amine", ("fluorinated_nitrogen_candidate", "etch_clean_candidate_review"), "small fluorinated amine candidate", None),
    SemiconductorAmineSeed("2,2-difluoroethylamine", "FC(F)CN", "primary", ("difluoroethyl",), "fluorinated_amine", ("fluorinated_nitrogen_candidate", "etch_clean_candidate_review"), "small fluorinated amine candidate", None),
    SemiconductorAmineSeed("2,2,2-trifluoroethylamine", "FC(F)(F)CN", "primary", ("trifluoroethyl",), "fluorinated_amine", ("fluorinated_nitrogen_candidate", "etch_clean_candidate_review"), "small trifluoroethyl amine candidate", "753-90-2"),
    SemiconductorAmineSeed("bis(2,2,2-trifluoroethyl)amine", "N(CC(F)(F)F)CC(F)(F)F", "secondary", ("trifluoroethyl", "trifluoroethyl"), "fluorinated_amine", ("fluorinated_nitrogen_candidate", "pfas_screening_candidate"), "fluorinated secondary amine candidate", None),
    SemiconductorAmineSeed("tris(2,2,2-trifluoroethyl)amine", "N(CC(F)(F)F)(CC(F)(F)F)CC(F)(F)F", "tertiary", ("trifluoroethyl", "trifluoroethyl", "trifluoroethyl"), "fluorinated_amine", ("fluorinated_nitrogen_candidate", "pfas_screening_candidate"), "fluorinated tertiary amine candidate", None),
    # Aminosilanes, silylamines, and amido precursor motifs.
    SemiconductorAmineSeed("dimethylaminosilane", "[SiH3]N(C)C", "amino_silane", ("dimethylamino",), "aminosilane_precursor", ("silicon_precursor_candidate", "ald_cvd_precursor_candidate"), "amino-silane precursor motif", None),
    SemiconductorAmineSeed("diethylaminosilane", "[SiH3]N(CC)CC", "amino_silane", ("diethylamino",), "aminosilane_precursor", ("silicon_precursor_candidate", "ald_cvd_precursor_candidate"), "amino-silane precursor motif", None),
    SemiconductorAmineSeed("bis(dimethylamino)silane", "[SiH2](N(C)C)N(C)C", "amino_silane", ("dimethylamino", "dimethylamino"), "aminosilane_precursor", ("silicon_precursor_candidate", "ald_cvd_precursor_candidate"), "bis(amino)silane precursor motif", None),
    SemiconductorAmineSeed("bis(diethylamino)silane", "[SiH2](N(CC)CC)N(CC)CC", "amino_silane", ("diethylamino", "diethylamino"), "aminosilane_precursor", ("silicon_precursor_candidate", "ald_cvd_precursor_candidate"), "bis(amino)silane precursor motif", None),
    SemiconductorAmineSeed("tris(dimethylamino)silane", "[SiH](N(C)C)(N(C)C)N(C)C", "amino_silane", ("dimethylamino", "dimethylamino", "dimethylamino"), "aminosilane_precursor", ("silicon_precursor_candidate", "ald_cvd_precursor_candidate"), "tris(amino)silane precursor motif", None),
    SemiconductorAmineSeed("tetrakis(dimethylamino)silane", "[Si](N(C)C)(N(C)C)(N(C)C)N(C)C", "amino_silane", ("dimethylamino", "dimethylamino", "dimethylamino", "dimethylamino"), "aminosilane_precursor", ("silicon_precursor_candidate", "ald_cvd_precursor_candidate"), "tetrakis(amino)silane precursor motif", None),
    SemiconductorAmineSeed("hexamethyldisilazane", "C[Si](C)(C)N[Si](C)(C)C", "silylamine", ("trimethylsilyl", "trimethylsilyl"), "silylamine_surface_modifier", ("surface_passivation_candidate", "silylation_candidate"), "silylamine surface chemistry motif", "999-97-3"),
    SemiconductorAmineSeed("tert-butylaminosilane", "[SiH3]NC(C)(C)C", "amino_silane", ("tert-butylamino",), "aminosilane_precursor", ("silicon_precursor_candidate", "sterically_hindered_ligand_candidate"), "bulky amino-silane precursor motif", None),
    SemiconductorAmineSeed("bis(tert-butylamino)silane", "[SiH2](NC(C)(C)C)NC(C)(C)C", "amino_silane", ("tert-butylamino", "tert-butylamino"), "aminosilane_precursor", ("silicon_precursor_candidate", "sterically_hindered_ligand_candidate"), "bulky bis(amino)silane precursor motif", None),
    SemiconductorAmineSeed("tris(dimethylamino)borane", "B(N(C)C)(N(C)C)N(C)C", "boron_amide", ("dimethylamino", "dimethylamino", "dimethylamino"), "boron_amide_precursor", ("boron_precursor_candidate", "ald_cvd_precursor_candidate"), "amino-boron precursor motif", None),
    SemiconductorAmineSeed("tetrakis(dimethylamido)titanium", "[Ti](N(C)C)(N(C)C)(N(C)C)N(C)C", "metal_amide", ("dimethylamido", "dimethylamido", "dimethylamido", "dimethylamido"), "metal_amide_precursor", ("titanium_precursor_candidate", "ald_cvd_precursor_candidate"), "metal amido precursor motif", "3275-24-9"),
]


class SemiconductorAmineGenerator:
    name = "SemiconductorAmineGenerator"

    def __init__(
        self,
        max_candidates: int = 80,
        max_heavy_atoms: int = 20,
        allowed_precursor_families: list[str] | None = None,
    ):
        self.max_candidates = int(max_candidates)
        self.max_heavy_atoms = int(max_heavy_atoms)
        self.allowed_precursor_families = {
            item.strip() for item in allowed_precursor_families or [] if item.strip()
        }

    def generate(self) -> Iterable[MoleculeCandidate]:
        count = 0
        seen: set[str] = set()
        for seed in SEMICONDUCTOR_AMINE_SEEDS:
            if count >= self.max_candidates:
                return
            if self.allowed_precursor_families and seed.precursor_family not in self.allowed_precursor_families:
                continue
            parsed = parse_smiles_lightweight(seed.smiles)
            if parsed is None or parsed.unsupported_tokens or parsed.heavy_atom_count > self.max_heavy_atoms:
                continue
            if seed.smiles in seen:
                continue
            seen.add(seed.smiles)
            count += 1
            yield self._candidate(seed, parsed.heavy_atom_count)

    def _candidate(self, seed: SemiconductorAmineSeed, heavy_atom_count: int) -> MoleculeCandidate:
        cid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"semiconductor_amine:{seed.name}:{seed.smiles}"))
        fluorinated_count = sum(1 for item in seed.substituents if "fluoro" in item)
        metadata = {
            "candidate_scope": "semiconductor_amine_curated",
            "amine_class": seed.amine_class,
            "amine_substituents": "; ".join(seed.substituents),
            "amine_substituent_count": str(len(seed.substituents)),
            "fluorinated_substituent_count": str(fluorinated_count),
            "contains_fluorinated_substituent": str(fluorinated_count > 0 or "F" in seed.smiles).lower(),
            "contains_unsaturated_substituent": "false",
            "contains_cyclic_substituent": str("1" in seed.smiles or seed.amine_class == "cyclic").lower(),
            "precursor_family": seed.precursor_family,
            "semiconductor_process_roles": "; ".join(seed.process_roles),
            "semiconductor_relevance_basis": seed.relevance_basis,
            "local_heavy_atom_count": str(heavy_atom_count),
        }
        return MoleculeCandidate(
            candidate_id=cid,
            source="semiconductor_amine",
            family="amine",
            input_name=seed.name,
            smiles=seed.smiles,
            cas=seed.cas,
            generation_rule="semiconductor_amine_curated",
            generation_score=0.88,
            metadata=metadata,
        )
