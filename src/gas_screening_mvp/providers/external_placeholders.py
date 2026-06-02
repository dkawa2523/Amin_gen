from __future__ import annotations

from typing import Iterable

from gas_screening_mvp.domain.models import ChemicalIdentity, PropertyCandidate, ReactionProbeSummary


class ChemeoProvider:
    name = "Chemeo"

    def __init__(self, enabled: bool = False):
        self.enabled = enabled

    def supports(self, property_name: str) -> bool:
        return self.enabled and property_name in {
            "normal_boiling_point",
            "normal_melting_point",
            "critical_temperature",
            "critical_pressure",
            "vapor_pressure",
        }

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        if not self.enabled:
            return []
        raise NotImplementedError("Add Chemeo credentials/terms-approved connector before enabling.")


class NistWebBookProvider:
    name = "NISTWebBook"

    def __init__(self, enabled: bool = False):
        self.enabled = enabled

    def supports(self, property_name: str) -> bool:
        return self.enabled and property_name in {
            "normal_boiling_point",
            "normal_melting_point",
            "vapor_pressure",
        }

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        if not self.enabled:
            return []
        raise NotImplementedError("Implement shortlist-only NIST Antoine extraction with usage review before enabling.")


class NistKineticsProvider:
    name = "NISTKinetics"

    def __init__(self, enabled: bool = False):
        self.enabled = enabled

    def probe(self, chemical: ChemicalIdentity, target_species: list[dict]) -> list[ReactionProbeSummary]:
        if not self.enabled:
            return []
        raise NotImplementedError("Implement shortlist-only NIST kinetics search/curation before enabling.")


class LxcatProvider:
    name = "LXCat"

    def __init__(self, enabled: bool = False):
        self.enabled = enabled

    def probe(self, chemical: ChemicalIdentity, target_species: list[dict]) -> list[ReactionProbeSummary]:
        if not self.enabled:
            return []
        raise NotImplementedError("Implement LXCat data export/import and EEDF/rate processing before enabling.")


class QdbProvider:
    name = "QDB"

    def __init__(self, enabled: bool = False):
        self.enabled = enabled

    def probe(self, chemical: ChemicalIdentity, target_species: list[dict]) -> list[ReactionProbeSummary]:
        if not self.enabled:
            return []
        raise NotImplementedError("Implement licensed QDB connector after contract/API review.")
