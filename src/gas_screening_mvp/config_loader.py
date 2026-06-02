from __future__ import annotations

from pathlib import Path
from copy import deepcopy
import yaml


DEFAULT_CONFIG = {
    "mode": "enrichment",
    "generation": {
        "max_total_candidates": 200,
        "enabled_generators": ["SeedListGenerator"],
        "include_formula_only": False,
        "amine_max_candidates": 500,
        "amine_max_carbons": 6,
        "amine_max_heavy_atoms": 12,
        "amine_allowed_substituents": [],
        "semiconductor_amine_max_candidates": 80,
        "semiconductor_amine_max_heavy_atoms": 20,
        "semiconductor_amine_allowed_precursor_families": [],
        "local_mutation_enabled": False,
        "local_mutation_max_candidates": 100,
        "pubchem_expansion_enabled": False,
        "pubchem_similarity_threshold": 90,
        "max_cids_per_seed": 25,
        "pubchem_expansion_rps": 1.0,
        "formula_search_enabled": False,
        "max_cids_per_formula": 10,
    },
    "prefilter": {
        "max_molecular_weight": 300,
        "max_heavy_atoms": 20,
        "allow_formula_only": False,
        "min_score": 0.0,
    },
    "providers": {
        "pubchem_enabled": False,
        "pugview_enabled": False,
        "cache_db": "./gas_screening_cache.sqlite",
        "curated_properties_csv": None,
        "gwp_csv": None,
        "pfas_list_csv": None,
        "reactivity_csv": None,
        "kinetics_csv": None,
        "max_api_candidates_per_run": 25,
    },
    "selection": {
        "required_properties": [
            "normal_melting_point",
            "normal_boiling_point",
            "critical_temperature",
            "critical_pressure",
            "gwp100_ar6",
        ],
        "target_vapor_temperatures_K": [298.15, 313.15, 333.15],
    },
    "fetch_planning": {
        "min_identity_api_score": 0.45,
        "min_enrichment_api_score": 0.85,
    },
    "supply_thresholds": {
        "high_vapor_pressure_kPa": 100.0,
        "bubbler_kPa": 10.0,
        "heated_source_kPa": 1.0,
    },
    "reaction_probe_targets": [
        {"species": "HF", "state": None},
        {"species": "O", "state": "O(3P)"},
        {"species": "O", "state": "O(1D)"},
        {"species": "F", "state": "atom"},
        {"species": "OH", "state": "radical"},
        {"species": "Cl", "state": "atom"},
        {"species": "e-", "state": "electron"},
    ],
}


def deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: str | Path | None = None) -> dict:
    cfg = deepcopy(DEFAULT_CONFIG)
    if path:
        with Path(path).open(encoding="utf-8") as f:
            user_cfg = yaml.safe_load(f) or {}
        user_cfg = normalize_config_aliases(user_cfg)
        cfg = deep_merge(cfg, user_cfg)
    return cfg


def normalize_config_aliases(user_cfg: dict) -> dict:
    """Accept the more descriptive exploration config shape used in notebooks/docs."""
    cfg = deepcopy(user_cfg or {})
    generation = cfg.get("generation")
    if not isinstance(generation, dict):
        return cfg

    if "mode" in generation and "mode" not in cfg:
        cfg["mode"] = generation["mode"]

    if "max_api_candidates_per_run" in generation:
        cfg.setdefault("providers", {})["max_api_candidates_per_run"] = generation["max_api_candidates_per_run"]

    amine = generation.get("amine")
    if isinstance(amine, dict):
        if "max_candidates" in amine:
            generation["amine_max_candidates"] = amine["max_candidates"]
        if "max_carbons" in amine:
            generation["amine_max_carbons"] = amine["max_carbons"]
        if "max_heavy_atoms" in amine:
            generation["amine_max_heavy_atoms"] = amine["max_heavy_atoms"]
        if "substituents" in amine:
            generation["amine_allowed_substituents"] = amine["substituents"]
        elif amine.get("include_fluorinated_substituents"):
            generation["amine_allowed_substituents"] = ["small_alkyl", "alkyl", "cyclic_alkyl", "fluorinated_alkyl"]

    fluorocarbon = generation.get("fluorocarbon")
    if isinstance(fluorocarbon, dict):
        if "include_formula_only" in fluorocarbon:
            generation["include_formula_only"] = fluorocarbon["include_formula_only"]
        if "carbon_range" in fluorocarbon:
            generation["fluorocarbon_carbon_range"] = fluorocarbon["carbon_range"]
        if "max_cids_per_formula" in fluorocarbon:
            generation["max_cids_per_formula"] = fluorocarbon["max_cids_per_formula"]

    semiconductor_amine = generation.get("semiconductor_amine")
    if isinstance(semiconductor_amine, dict):
        if "max_candidates" in semiconductor_amine:
            generation["semiconductor_amine_max_candidates"] = semiconductor_amine["max_candidates"]
        if "max_heavy_atoms" in semiconductor_amine:
            generation["semiconductor_amine_max_heavy_atoms"] = semiconductor_amine["max_heavy_atoms"]
        if "allowed_precursor_families" in semiconductor_amine:
            generation["semiconductor_amine_allowed_precursor_families"] = semiconductor_amine["allowed_precursor_families"]

    pubchem_expansion = generation.get("pubchem_expansion")
    if isinstance(pubchem_expansion, dict):
        if "enabled" in pubchem_expansion:
            generation["pubchem_expansion_enabled"] = pubchem_expansion["enabled"]
        if "threshold" in pubchem_expansion:
            generation["pubchem_similarity_threshold"] = pubchem_expansion["threshold"]
        if "max_cids_per_seed" in pubchem_expansion:
            generation["max_cids_per_seed"] = pubchem_expansion["max_cids_per_seed"]
        if "rps" in pubchem_expansion:
            generation["pubchem_expansion_rps"] = pubchem_expansion["rps"]

    local_mutation = generation.get("local_mutation")
    if isinstance(local_mutation, dict):
        if "enabled" in local_mutation:
            generation["local_mutation_enabled"] = local_mutation["enabled"]
        if "max_candidates_per_seed" in local_mutation:
            generation["local_mutation_max_candidates"] = local_mutation["max_candidates_per_seed"]

    if "planning" in cfg and "fetch_planning" not in cfg:
        cfg["fetch_planning"] = cfg["planning"]

    return cfg
