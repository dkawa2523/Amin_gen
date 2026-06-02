from __future__ import annotations

import json
import re
from typing import Iterable

import requests

from gas_screening_mvp.domain.api import ApiRequest
from gas_screening_mvp.domain.models import ChemicalIdentity, ClassificationRecord, PropertyCandidate
from gas_screening_mvp.storage.cache import SqliteApiCache
from gas_screening_mvp.planning.rate_limiter import SyncRateLimiter
from gas_screening_mvp.providers.reactivity import PHYSICAL_H_CODES


class PubChemPugViewProvider:
    """Shortlist-only PubChem PUG-View provider for GHS H-codes and text fallback.

    This is intentionally limited and conservative. It should not be used as a
    high-volume property source.
    """

    name = "PubChemPugView"
    BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound"

    def __init__(self, cache: SqliteApiCache, enabled: bool = False, rps: float = 1.0, timeout: int = 30):
        self.cache = cache
        self.enabled = enabled
        self.timeout = timeout
        self.limiter = SyncRateLimiter({self.name: rps})

    def classify(self, chemical: ChemicalIdentity) -> list[ClassificationRecord]:
        if not self.enabled or not chemical.pubchem_cid:
            return []
        data = self._fetch_heading(chemical.pubchem_cid, "GHS Classification")
        if not data:
            return []
        text = json.dumps(data, ensure_ascii=False)
        hcodes = sorted(set(re.findall(r"\bH\d{3}\b", text)) & PHYSICAL_H_CODES)
        return [
            ClassificationRecord(
                chemical.chemical_id,
                "ghs_physical_h_code",
                h,
                "PubChem PUG-View GHS Classification text match",
                self.name,
                confidence="C",
            )
            for h in hcodes
        ]

    def supports(self, property_name: str) -> bool:
        return property_name in {"normal_boiling_point", "normal_melting_point", "vapor_pressure"}

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        # MVP: leave experimental-text property extraction as a controlled Phase 1B task.
        # Keeping this empty prevents accidental use of weak text-extracted values.
        return []

    def _fetch_heading(self, cid: int, heading: str) -> dict | None:
        endpoint = f"{self.BASE}/{cid}/JSON"
        params = {"heading": heading}
        req = ApiRequest(self.name, endpoint, "GET", params, None)
        cached = self.cache.get(req.signature())
        if cached:
            text, status = cached["response_body"], int(cached["status_code"])
        else:
            self.limiter.wait(self.name)
            res = requests.get(endpoint, params=params, timeout=self.timeout)
            self.cache.put(req, res.status_code, res.text, ttl_days=90)
            text, status = res.text, res.status_code
        if status != 200:
            return None
        try:
            return json.loads(text)
        except Exception:
            return None
