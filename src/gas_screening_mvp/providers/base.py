from __future__ import annotations

from typing import Protocol, Iterable

from gas_screening_mvp.domain.models import (
    NormalizedMolecule,
    ChemicalIdentity,
    PropertyCandidate,
    ClassificationRecord,
    ReactionProbeSummary,
)


class IdentityProvider(Protocol):
    name: str

    def resolve(self, molecule: NormalizedMolecule) -> list[ChemicalIdentity]:
        ...


class PropertyProvider(Protocol):
    name: str

    def supports(self, property_name: str) -> bool:
        ...

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        ...


class ClassificationProvider(Protocol):
    name: str

    def classify(self, chemical: ChemicalIdentity) -> list[ClassificationRecord]:
        ...


class KineticsProvider(Protocol):
    name: str

    def probe(self, chemical: ChemicalIdentity, target_species: list[dict]) -> list[ReactionProbeSummary]:
        ...
