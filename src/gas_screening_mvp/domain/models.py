"""Domain models for the semiconductor gas screening MVP.

The core design principle is to keep every fetched or calculated value as a
PropertyCandidate, then separately store the SelectedProperty that appears in
engineering-facing outputs.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Literal


def utcnow() -> datetime:
    return datetime.now(timezone.utc)


CandidateSource = Literal[
    "seed",
    "pubchem_similarity",
    "pubchem_substructure",
    "pubchem_formula",
    "template_generation",
    "semiconductor_amine",
    "formula_generation",
    "local_mutation",
    "manual",
]


@dataclass(frozen=True)
class MoleculeCandidate:
    candidate_id: str
    source: CandidateSource
    family: str
    input_name: str | None = None
    smiles: str | None = None
    formula: str | None = None
    cas: str | None = None
    generation_rule: str = "manual"
    parent_candidate_id: str | None = None
    generation_score: float = 0.0
    metadata: dict[str, str] = field(default_factory=dict)


StructureStatus = Literal[
    "valid",
    "invalid",
    "salt",
    "mixture",
    "radical",
    "unsupported",
    "manual_review_required",
]


@dataclass(frozen=True)
class NormalizedMolecule:
    candidate_id: str
    source: str
    family: str
    input_name: str | None
    cas: str | None
    canonical_smiles: str | None
    isomeric_smiles: str | None
    standard_inchi: str | None
    standard_inchikey: str | None
    formula: str | None
    molecular_weight: float | None
    heavy_atom_count: int | None
    element_symbols: tuple[str, ...]
    structure_status: StructureStatus
    status_reason: str | None = None
    metadata: dict[str, str] = field(default_factory=dict)


IdentityStatus = Literal[
    "resolved",
    "ambiguous",
    "unresolved",
    "manual_review_required",
]


@dataclass(frozen=True)
class ChemicalIdentity:
    chemical_id: str
    candidate_id: str
    preferred_name: str | None
    cas: str | None
    pubchem_cid: int | None
    formula: str | None
    molecular_weight: float | None
    canonical_smiles: str | None
    isomeric_smiles: str | None
    inchi: str | None
    inchikey: str | None
    identity_status: IdentityStatus
    confidence: float
    source: str


@dataclass(frozen=True)
class PropertyCandidate:
    chemical_id: str
    property_name: str
    value_num: float | None = None
    value_text: str | None = None
    unit: str | None = None
    temperature_K: float | None = None
    pressure_Pa: float | None = None
    phase: str | None = None
    source: str = ""
    method: str | None = None
    source_version: str | None = None
    reference: str | None = None
    valid_temperature_min_K: float | None = None
    valid_temperature_max_K: float | None = None
    is_estimated: bool = False
    retrieved_at: datetime = field(default_factory=utcnow)
    raw_id: str | None = None
    quality_hint: str | None = None


SelectedStatus = Literal[
    "selected",
    "missing",
    "not_applicable",
    "outside_range",
    "conflict",
    "manual_review_required",
]


@dataclass(frozen=True)
class SelectedProperty:
    chemical_id: str
    property_name: str
    value_num: float | None
    value_text: str | None
    unit: str | None
    status: SelectedStatus
    quality_rank: Literal["A", "B", "C", "D", "Missing", "N/A", "Conflict"]
    selected_source: str | None
    selection_reason: str
    selected_at: datetime = field(default_factory=utcnow)


@dataclass(frozen=True)
class ClassificationRecord:
    chemical_id: str
    classification_type: Literal[
        "pfas_flag",
        "pfas_basis",
        "pfas_definition",
        "pfas_motif",
        "pfas_list_hit",
        "persistence_screen",
        "persistence_basis",
        "reactive_group",
        "reactivity_flag",
        "ghs_physical_h_code",
    ]
    value: str
    basis: str
    source: str
    source_version: str | None = None
    confidence: Literal["A", "B", "C", "D"] = "C"
    retrieved_at: datetime = field(default_factory=utcnow)


@dataclass(frozen=True)
class ReactionProbeSummary:
    chemical_id: str
    target_species: str
    target_species_state: str | None
    phase: Literal["thermal_gas", "plasma", "surface", "unknown"]
    availability: Literal["available", "partial", "none", "not_checked"]
    sources: tuple[str, ...] = ()
    comment: str | None = None


@dataclass(frozen=True)
class SummaryRow:
    values: dict[str, object]
