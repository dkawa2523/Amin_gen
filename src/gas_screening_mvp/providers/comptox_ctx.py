from __future__ import annotations

from gas_screening_mvp.domain.models import ChemicalIdentity, ClassificationRecord


class CompToxCtxProvider:
    """Placeholder for EPA CTX/CompTox API enrichment.

    MVP production recommendation: prefer periodic download/local mirroring of
    PFAS/list memberships where possible, then use CTX API only for unresolved
    or high-priority molecules. This class is intentionally not wired into the
    default pipeline until endpoint contracts and API keys are configured by the
    deployment owner.
    """

    name = "CompToxCTX"

    def __init__(self, api_key: str | None = None, enabled: bool = False):
        self.api_key = api_key
        self.enabled = enabled

    def classify(self, chemical: ChemicalIdentity) -> list[ClassificationRecord]:
        if not self.enabled:
            return []
        raise NotImplementedError("Configure CTX endpoint URLs, API key handling, and response mapping before enabling.")
