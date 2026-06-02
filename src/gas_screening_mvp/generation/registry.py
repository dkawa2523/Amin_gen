from __future__ import annotations

from pathlib import Path
from typing import Iterable

from gas_screening_mvp.domain.models import MoleculeCandidate
from gas_screening_mvp.generation.seed_generator import SeedListGenerator
from gas_screening_mvp.generation.amine_template import AmineTemplateGenerator
from gas_screening_mvp.generation.fluorocarbon_formula import FluorocarbonFormulaGenerator
from gas_screening_mvp.generation.local_mutation import LocalMutationGenerator
from gas_screening_mvp.generation.pubchem_expansion import PubChemExpansionGenerator
from gas_screening_mvp.generation.semiconductor_amine import SemiconductorAmineGenerator
from gas_screening_mvp.storage.cache import SqliteApiCache


def load_generators(
    config: dict,
    input_csv: str | Path | None = None,
    cache: SqliteApiCache | None = None,
    mode: str = "enrichment",
    remote_enabled: bool = False,
    dry_run: bool = False,
):
    enabled = set(config.get("enabled_generators", ["SeedListGenerator"]))
    gens = []
    seed_candidates: list[MoleculeCandidate] = []
    if "SeedListGenerator" in enabled:
        seed_gen = SeedListGenerator(input_csv)
        seed_candidates = list(seed_gen.generate())
        gens.append(seed_candidates)
    if "SemiconductorAmineGenerator" in enabled:
        gens.append(SemiconductorAmineGenerator(
            max_candidates=config.get("semiconductor_amine_max_candidates", 80),
            max_heavy_atoms=config.get("semiconductor_amine_max_heavy_atoms", config.get("amine_max_heavy_atoms", 20)),
            allowed_precursor_families=config.get("semiconductor_amine_allowed_precursor_families") or None,
        ))
    if "AmineTemplateGenerator" in enabled:
        gens.append(AmineTemplateGenerator(
            max_candidates=config.get("amine_max_candidates", 500),
            max_carbons=config.get("amine_max_carbons", 6),
            max_heavy_atoms=config.get("amine_max_heavy_atoms", 12),
            allowed_substituents=config.get("amine_allowed_substituents") or None,
        ))
    if "FluorocarbonFormulaGenerator" in enabled:
        include_formula_only = bool(config.get("include_formula_only", False)) and mode == "exploration"
        gens.append(FluorocarbonFormulaGenerator(
            include_formula_only=include_formula_only,
            carbon_range=config.get("fluorocarbon_carbon_range"),
        ))
    if "LocalMutationGenerator" in enabled and mode == "exploration":
        gens.append(LocalMutationGenerator(
            seeds=seed_candidates,
            enabled=bool(config.get("local_mutation_enabled", False)),
            max_candidates=config.get("local_mutation_max_candidates", 100),
            max_heavy_atoms=config.get("local_mutation_max_heavy_atoms", config.get("amine_max_heavy_atoms", 20)),
            allowed_families=config.get("local_mutation_allowed_families") or None,
        ))
    if (
        "PubChemExpansionGenerator" in enabled
        and mode == "exploration"
        and remote_enabled
        and not dry_run
        and cache is not None
        and bool(config.get("pubchem_expansion_enabled", False))
    ):
        seed_smiles = [c.smiles for c in seed_candidates if c.smiles]
        gens.append(PubChemExpansionGenerator(
            seed_smiles=seed_smiles,
            cache=cache,
            enabled=True,
            threshold=int(config.get("pubchem_similarity_threshold", 90)),
            max_cids_per_seed=int(config.get("max_cids_per_seed", 25)),
            rps=float(config.get("pubchem_expansion_rps", 1.0)),
        ))
    return gens


def generate_candidates(
    config: dict,
    input_csv: str | Path | None = None,
    cache: SqliteApiCache | None = None,
    mode: str = "enrichment",
    remote_enabled: bool = False,
    dry_run: bool = False,
) -> list[MoleculeCandidate]:
    out: list[MoleculeCandidate] = []
    max_total = int(config.get("max_total_candidates", 50000))
    for gen in load_generators(config, input_csv, cache=cache, mode=mode, remote_enabled=remote_enabled, dry_run=dry_run):
        source = gen if isinstance(gen, list) else gen.generate()
        for cand in source:
            out.append(cand)
            if len(out) >= max_total:
                return out
    return out
