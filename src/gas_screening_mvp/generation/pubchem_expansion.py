from __future__ import annotations

import uuid
from typing import Iterable

import requests

from gas_screening_mvp.domain.models import MoleculeCandidate
from gas_screening_mvp.domain.api import ApiRequest
from gas_screening_mvp.storage.cache import SqliteApiCache
from gas_screening_mvp.planning.rate_limiter import SyncRateLimiter


class PubChemExpansionGenerator:
    """Optional generator that expands seed SMILES through PubChem.

    The generator is not enabled by default. It should be used in exploration
    mode with conservative limits. Returned candidates initially contain CID as
    metadata; structural properties should be fetched later by PUG-REST batch.
    """

    name = "PubChemExpansionGenerator"
    BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(
        self,
        seed_smiles: list[str],
        cache: SqliteApiCache,
        enabled: bool = False,
        threshold: int = 90,
        max_cids_per_seed: int = 50,
        rps: float = 1.0,
    ):
        self.seed_smiles = seed_smiles
        self.cache = cache
        self.enabled = enabled
        self.threshold = threshold
        self.max_cids_per_seed = max_cids_per_seed
        self.limiter = SyncRateLimiter({self.name: rps})

    def generate(self) -> Iterable[MoleculeCandidate]:
        if not self.enabled:
            return []
        out: list[MoleculeCandidate] = []
        for smi in self.seed_smiles:
            cids = self._similarity_cids(smi)[: self.max_cids_per_seed]
            for cid in cids:
                candidate_id = str(uuid.uuid5(uuid.NAMESPACE_URL, f"pubchem_similarity:{smi}:{cid}"))
                out.append(MoleculeCandidate(
                    candidate_id=candidate_id,
                    source="pubchem_similarity",
                    family="pubchem_expansion",
                    input_name=f"PubChem CID {cid}",
                    smiles=None,
                    formula=None,
                    generation_rule=f"fastsimilarity_2d_threshold_{self.threshold}",
                    parent_candidate_id=None,
                    generation_score=0.75,
                    metadata={"pubchem_cid": str(cid), "seed_smiles": smi},
                ))
        return out

    def _similarity_cids(self, smiles: str) -> list[int]:
        if self.cache.has_negative(self.name, "similarity_smiles", smiles):
            return []
        endpoint = f"{self.BASE}/compound/fastsimilarity_2d/smiles/cids/txt"
        req = ApiRequest(self.name, endpoint, "POST", {"Threshold": self.threshold}, {"smiles": smiles})
        cached = self.cache.get(req.signature())
        if cached:
            text = cached["response_body"]
            status = int(cached["status_code"])
        else:
            self.limiter.wait(self.name)
            res = requests.post(endpoint, params={"Threshold": self.threshold}, data={"smiles": smiles}, timeout=60)
            self.cache.put(req, res.status_code, res.text, ttl_days=30)
            text, status = res.text, res.status_code
        if status != 200:
            self.cache.put_negative(self.name, "similarity_smiles", smiles, f"status_{status}", ttl_days=30)
            return []
        cids = []
        for token in text.split():
            try:
                cids.append(int(token))
            except ValueError:
                pass
        if not cids:
            self.cache.put_negative(self.name, "similarity_smiles", smiles, "no_cids", ttl_days=30)
        return cids
