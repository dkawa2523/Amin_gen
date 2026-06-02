from __future__ import annotations

import csv
import io
import json
import uuid
import re
from typing import Iterable

import requests

from gas_screening_mvp.domain.api import ApiRequest
from gas_screening_mvp.domain.models import NormalizedMolecule, ChemicalIdentity, PropertyCandidate
from gas_screening_mvp.storage.cache import SqliteApiCache
from gas_screening_mvp.planning.rate_limiter import SyncRateLimiter


class PubChemPugRestProvider:
    """PubChem PUG-REST identity and basic property resolver.

    Uses cache + rate limit. Network access is controlled by `enabled` and
    should remain disabled in unit tests/offline runs.
    """

    name = "PubChemPugRest"
    BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    BASIC_PROPERTIES = [
        "MolecularFormula",
        "MolecularWeight",
        "CanonicalSMILES",
        "IsomericSMILES",
        "InChI",
        "InChIKey",
        "IUPACName",
    ]

    def __init__(self, cache: SqliteApiCache, enabled: bool = False, rps: float = 3.0, timeout: int = 30):
        self.cache = cache
        self.enabled = enabled
        self.timeout = timeout
        self.limiter = SyncRateLimiter({self.name: rps})

    def resolve(self, molecule: NormalizedMolecule) -> list[ChemicalIdentity]:
        if not self.enabled:
            return []
        cid = None
        if molecule.metadata.get("pubchem_cid"):
            try:
                cid = int(molecule.metadata["pubchem_cid"])
            except ValueError:
                cid = None
        if cid is None and molecule.standard_inchikey:
            cid = self.cid_from_inchikey(molecule.standard_inchikey)
        if cid is None and molecule.canonical_smiles:
            cid = self.cid_from_smiles(molecule.canonical_smiles)
        if cid is None and molecule.input_name:
            cid = self.cid_from_name(molecule.input_name)
        if cid is None:
            return []
        rows = self.properties_for_cids([cid])
        if not rows:
            return []
        row = rows[0]
        key = row.get("InChIKey") or molecule.standard_inchikey or str(cid)
        return [
            ChemicalIdentity(
                chemical_id=str(uuid.uuid5(uuid.NAMESPACE_URL, f"pubchem:{key}")),
                candidate_id=molecule.candidate_id,
                preferred_name=row.get("IUPACName") or molecule.input_name,
                cas=molecule.cas,
                pubchem_cid=int(row["CID"]),
                formula=row.get("MolecularFormula") or molecule.formula,
                molecular_weight=_float_or_none(row.get("MolecularWeight")) or molecule.molecular_weight,
                canonical_smiles=row.get("CanonicalSMILES") or molecule.canonical_smiles,
                isomeric_smiles=row.get("IsomericSMILES") or molecule.isomeric_smiles,
                inchi=row.get("InChI") or molecule.standard_inchi,
                inchikey=row.get("InChIKey") or molecule.standard_inchikey,
                identity_status="resolved",
                confidence=0.9,
                source=self.name,
            )
        ]

    def resolve_formula(
        self,
        molecule: NormalizedMolecule,
        max_cids_per_formula: int = 10,
        allowed_elements: set[str] | None = None,
        max_heavy_atoms: int | None = None,
        max_molecular_weight: float | None = None,
    ) -> list[ChemicalIdentity]:
        if not self.enabled or not molecule.formula:
            return []
        formula_counts = _parse_formula(molecule.formula)
        if not formula_counts:
            return []
        if allowed_elements and not set(formula_counts).issubset(allowed_elements):
            return []
        if max_heavy_atoms is not None and _heavy_atom_count(formula_counts) > max_heavy_atoms:
            return []
        cids = self.cids_from_formula(molecule.formula, max_cids=max_cids_per_formula)
        if not cids:
            return []
        rows = self.properties_for_cids(cids)
        hits: list[ChemicalIdentity] = []
        for row in rows:
            row_formula = row.get("MolecularFormula") or ""
            row_counts = _parse_formula(row_formula)
            if row_counts != formula_counts:
                continue
            if allowed_elements and not set(row_counts).issubset(allowed_elements):
                continue
            if max_heavy_atoms is not None and _heavy_atom_count(row_counts) > max_heavy_atoms:
                continue
            mw = _float_or_none(row.get("MolecularWeight"))
            if max_molecular_weight is not None and mw is not None and mw > max_molecular_weight:
                continue
            cid = int(row["CID"])
            key = row.get("InChIKey") or str(cid)
            hits.append(ChemicalIdentity(
                chemical_id=str(uuid.uuid5(uuid.NAMESPACE_URL, f"pubchem:{key}")),
                candidate_id=molecule.candidate_id,
                preferred_name=row.get("IUPACName") or molecule.input_name,
                cas=molecule.cas,
                pubchem_cid=cid,
                formula=row_formula or molecule.formula,
                molecular_weight=mw or molecule.molecular_weight,
                canonical_smiles=row.get("CanonicalSMILES") or molecule.canonical_smiles,
                isomeric_smiles=row.get("IsomericSMILES") or molecule.isomeric_smiles,
                inchi=row.get("InChI") or molecule.standard_inchi,
                inchikey=row.get("InChIKey") or molecule.standard_inchikey,
                identity_status="resolved",
                confidence=0.82,
                source=self.name,
            ))
        if not hits:
            self.cache.put_negative(self.name, "formula", molecule.formula, "no_filtered_hits", ttl_days=30)
        return hits

    def cid_from_inchikey(self, inchikey: str) -> int | None:
        if self.cache.has_negative(self.name, "inchikey", inchikey):
            return None
        endpoint = f"{self.BASE}/compound/inchikey/{inchikey}/cids/txt"
        req = ApiRequest(self.name, endpoint, "GET", {}, None)
        text, status = self._cached_request(req)
        if status != 200:
            self.cache.put_negative(self.name, "inchikey", inchikey, f"status_{status}", ttl_days=90)
            return None
        return _parse_first_int(text)

    def cid_from_name(self, name: str) -> int | None:
        if self.cache.has_negative(self.name, "name", name):
            return None
        endpoint = f"{self.BASE}/compound/name/{requests.utils.quote(name)}/cids/txt"
        req = ApiRequest(self.name, endpoint, "GET", {}, None)
        text, status = self._cached_request(req)
        if status != 200:
            self.cache.put_negative(self.name, "name", name, f"status_{status}", ttl_days=90)
            return None
        return _parse_first_int(text)

    def cid_from_smiles(self, smiles: str) -> int | None:
        if self.cache.has_negative(self.name, "smiles", smiles):
            return None
        endpoint = f"{self.BASE}/compound/smiles/cids/txt"
        req = ApiRequest(self.name, endpoint, "POST", {}, {"smiles": smiles})
        text, status = self._cached_request(req)
        if status != 200:
            self.cache.put_negative(self.name, "smiles", smiles, f"status_{status}", ttl_days=90)
            return None
        return _parse_first_int(text)

    def cids_from_formula(self, formula: str, max_cids: int = 10) -> list[int]:
        if self.cache.has_negative(self.name, "formula", formula):
            return []
        endpoint = f"{self.BASE}/compound/fastformula/{requests.utils.quote(formula)}/cids/txt"
        req = ApiRequest(self.name, endpoint, "GET", {"MaxRecords": max_cids}, None)
        text, status = self._cached_request(req)
        if status != 200:
            self.cache.put_negative(self.name, "formula", formula, f"status_{status}", ttl_days=30)
            return []
        cids: list[int] = []
        for token in str(text).split():
            try:
                cids.append(int(token))
            except ValueError:
                continue
            if len(cids) >= max_cids:
                break
        if not cids:
            self.cache.put_negative(self.name, "formula", formula, "no_cids", ttl_days=30)
        return cids

    def properties_for_cids(self, cids: list[int]) -> list[dict[str, str]]:
        if not cids:
            return []
        cid_text = ",".join(map(str, cids))
        prop_text = ",".join(self.BASIC_PROPERTIES)
        endpoint = f"{self.BASE}/compound/cid/{cid_text}/property/{prop_text}/csv"
        req = ApiRequest(self.name, endpoint, "GET", {}, None)
        text, status = self._cached_request(req)
        if status != 200:
            return []
        return list(csv.DictReader(io.StringIO(text)))

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        # PUG-REST property endpoint is not used for thermal properties here.
        return []

    def supports(self, property_name: str) -> bool:
        return False

    def _cached_request(self, req: ApiRequest) -> tuple[str, int]:
        cached = self.cache.get(req.signature())
        if cached:
            return cached["response_body"], int(cached["status_code"])

        # Negative cache check uses the body or endpoint last component only in caller.
        self.limiter.wait(self.name)
        if req.method == "GET":
            res = requests.get(req.endpoint, params=req.params or None, timeout=self.timeout)
        else:
            res = requests.post(req.endpoint, params=req.params or None, data=req.body or {}, timeout=self.timeout)
        self.cache.put(req, res.status_code, res.text, ttl_days=180)
        return res.text, res.status_code


def _parse_first_int(text: str | None) -> int | None:
    if not text:
        return None
    for token in str(text).strip().split():
        try:
            return int(token)
        except ValueError:
            continue
    return None


def _float_or_none(value: object) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except Exception:
        return None


def _parse_formula(formula: str | None) -> dict[str, int]:
    if not formula:
        return {}
    counts: dict[str, int] = {}
    pos = 0
    for match in re.finditer(r"([A-Z][a-z]?)(\d*)", formula.replace(" ", "")):
        if match.start() != pos:
            return {}
        symbol = match.group(1)
        count = int(match.group(2) or "1")
        counts[symbol] = counts.get(symbol, 0) + count
        pos = match.end()
    return counts if pos == len(formula.replace(" ", "")) else {}


def _heavy_atom_count(counts: dict[str, int]) -> int:
    return sum(count for symbol, count in counts.items() if symbol != "H")
