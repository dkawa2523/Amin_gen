from __future__ import annotations

import json
import uuid

import requests

from gas_screening_mvp.domain.api import ApiRequest
from gas_screening_mvp.domain.models import NormalizedMolecule, ChemicalIdentity
from gas_screening_mvp.storage.cache import SqliteApiCache
from gas_screening_mvp.planning.rate_limiter import SyncRateLimiter


class OpsinIdentityProvider:
    """Optional OPSIN resolver for systematic organic names.

    Disabled by default. Use only when PubChem/name resolution fails and the
    input looks like a systematic organic name.
    """

    name = "OPSIN"
    BASE = "https://opsin.ch.cam.ac.uk/opsin"

    def __init__(self, cache: SqliteApiCache, enabled: bool = False, rps: float = 0.5, timeout: int = 30):
        self.cache = cache
        self.enabled = enabled
        self.timeout = timeout
        self.limiter = SyncRateLimiter({self.name: rps})

    def resolve(self, molecule: NormalizedMolecule) -> list[ChemicalIdentity]:
        if not self.enabled or not molecule.input_name:
            return []
        endpoint = f"{self.BASE}/{requests.utils.quote(molecule.input_name)}.json"
        req = ApiRequest(self.name, endpoint, "GET", {}, None)
        cached = self.cache.get(req.signature())
        if cached:
            text, status = cached["response_body"], int(cached["status_code"])
        else:
            self.limiter.wait(self.name)
            res = requests.get(endpoint, timeout=self.timeout)
            self.cache.put(req, res.status_code, res.text, ttl_days=180)
            text, status = res.text, res.status_code
        if status != 200:
            return []
        try:
            data = json.loads(text)
        except Exception:
            return []
        if data.get("status") != "SUCCESS":
            return []
        inchikey = data.get("stdInChIKey")
        smiles = data.get("smiles")
        inchi = data.get("stdInChI")
        key = inchikey or smiles or molecule.input_name
        return [ChemicalIdentity(
            chemical_id=str(uuid.uuid5(uuid.NAMESPACE_URL, f"opsin:{key}")),
            candidate_id=molecule.candidate_id,
            preferred_name=molecule.input_name,
            cas=molecule.cas,
            pubchem_cid=None,
            formula=molecule.formula,
            molecular_weight=molecule.molecular_weight,
            canonical_smiles=smiles or molecule.canonical_smiles,
            isomeric_smiles=smiles or molecule.isomeric_smiles,
            inchi=inchi or molecule.standard_inchi,
            inchikey=inchikey or molecule.standard_inchikey,
            identity_status="resolved" if inchikey else "manual_review_required",
            confidence=0.7 if inchikey else 0.4,
            source=self.name,
        )]
