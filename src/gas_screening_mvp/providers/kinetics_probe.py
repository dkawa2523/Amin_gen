from __future__ import annotations

import csv
from pathlib import Path

from gas_screening_mvp.domain.models import ChemicalIdentity, ReactionProbeSummary


class LocalKineticsProbeProvider:
    """Coverage-only kinetics probe provider.

    CSV columns: inchikey,cas,name,target_species,target_species_state,phase,
    availability,sources,comment
    """

    name = "LocalKineticsProbe"

    def __init__(self, csv_path: str | Path | None = None):
        self.rows: list[dict] = []
        if csv_path and Path(csv_path).exists():
            with Path(csv_path).open(newline="", encoding="utf-8-sig") as f:
                self.rows = list(csv.DictReader(f))

    def probe(self, chemical: ChemicalIdentity, target_species: list[dict]) -> list[ReactionProbeSummary]:
        out: list[ReactionProbeSummary] = []
        for target in target_species:
            species = target.get("species")
            state = target.get("state")
            matches = [r for r in self.rows if self._matches(r, chemical, species, state)]
            if matches:
                best = matches[0]
                out.append(ReactionProbeSummary(
                    chemical_id=chemical.chemical_id,
                    target_species=species,
                    target_species_state=state,
                    phase=best.get("phase") or "unknown",
                    availability=best.get("availability") or "partial",
                    sources=tuple(_split(best.get("sources"))),
                    comment=best.get("comment") or None,
                ))
            else:
                out.append(ReactionProbeSummary(
                    chemical_id=chemical.chemical_id,
                    target_species=species,
                    target_species_state=state,
                    phase="unknown",
                    availability="not_checked",
                    sources=(),
                    comment="No local kinetics coverage table match; NIST/LXCat/QDB connectors are Phase 2/shortlist tasks.",
                ))
        return out

    def _matches(self, row: dict, chemical: ChemicalIdentity, species: str, state: str | None) -> bool:
        if row.get("target_species") != species:
            return False
        if (row.get("target_species_state") or None) != state:
            return False
        if row.get("inchikey") and chemical.inchikey and row["inchikey"] == chemical.inchikey:
            return True
        if row.get("cas") and chemical.cas and row["cas"] == chemical.cas:
            return True
        if row.get("name") and chemical.preferred_name and row["name"].lower() == chemical.preferred_name.lower():
            return True
        return False


def _split(text: str | None) -> list[str]:
    if not text:
        return []
    return [x.strip() for x in text.replace(";", ",").split(",") if x.strip()]
