from __future__ import annotations

import sqlite3
from pathlib import Path
from datetime import datetime, timedelta, timezone
import json
import hashlib

from gas_screening_mvp.domain.api import ApiRequest


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _dt_to_str(dt: datetime | None) -> str | None:
    return dt.isoformat() if dt else None


class SqliteApiCache:
    def __init__(self, db_path: str | Path):
        self.db_path = Path(db_path)
        self.stats = {
            "cache_hits": 0,
            "negative_cache_hits": 0,
            "pubchem_requests": 0,
            "pugview_requests": 0,
        }
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init()

    def _connect(self):
        return sqlite3.connect(self.db_path)

    def _init(self) -> None:
        with self._connect() as con:
            con.execute(
                """
                CREATE TABLE IF NOT EXISTS api_request_cache (
                    signature TEXT PRIMARY KEY,
                    provider TEXT NOT NULL,
                    endpoint TEXT NOT NULL,
                    method TEXT NOT NULL,
                    params_json TEXT,
                    body_json TEXT,
                    status_code INTEGER,
                    response_body TEXT,
                    retrieved_at TEXT,
                    expires_at TEXT,
                    error_type TEXT
                )
                """
            )
            con.execute(
                """
                CREATE TABLE IF NOT EXISTS negative_cache (
                    provider TEXT NOT NULL,
                    query_type TEXT NOT NULL,
                    query_hash TEXT NOT NULL,
                    query_value TEXT,
                    negative_reason TEXT,
                    created_at TEXT,
                    expires_at TEXT,
                    PRIMARY KEY(provider, query_type, query_hash)
                )
                """
            )
            con.commit()

    def get(self, signature: str) -> dict | None:
        with self._connect() as con:
            cur = con.execute(
                "SELECT signature, provider, status_code, response_body, expires_at FROM api_request_cache WHERE signature=?",
                (signature,),
            )
            row = cur.fetchone()
            if not row:
                return None
            expires_at = datetime.fromisoformat(row[4]) if row[4] else None
            if expires_at and expires_at < _now():
                return None
            self.stats["cache_hits"] += 1
            return {
                "signature": row[0],
                "provider": row[1],
                "status_code": row[2],
                "response_body": row[3],
                "expires_at": row[4],
            }

    def put(self, req: ApiRequest, status_code: int, response_body: str, ttl_days: int | None = None, error_type: str | None = None) -> None:
        if req.provider == "PubChemPugView":
            self.stats["pugview_requests"] += 1
        elif req.provider.startswith("PubChem"):
            self.stats["pubchem_requests"] += 1
        now = _now()
        expires = now + timedelta(days=ttl_days) if ttl_days else None
        with self._connect() as con:
            con.execute(
                """
                INSERT OR REPLACE INTO api_request_cache
                (signature, provider, endpoint, method, params_json, body_json, status_code, response_body, retrieved_at, expires_at, error_type)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    req.signature(),
                    req.provider,
                    req.endpoint,
                    req.method,
                    json.dumps(req.params, sort_keys=True, ensure_ascii=False),
                    json.dumps(req.body, sort_keys=True, ensure_ascii=False) if req.body is not None else None,
                    status_code,
                    response_body,
                    _dt_to_str(now),
                    _dt_to_str(expires),
                    error_type,
                ),
            )
            con.commit()

    def put_negative(self, provider: str, query_type: str, query_value: str, reason: str, ttl_days: int = 90) -> None:
        now = _now()
        expires = now + timedelta(days=ttl_days)
        qhash = hashlib.sha256(query_value.encode("utf-8")).hexdigest()
        with self._connect() as con:
            con.execute(
                """
                INSERT OR REPLACE INTO negative_cache
                (provider, query_type, query_hash, query_value, negative_reason, created_at, expires_at)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (provider, query_type, qhash, query_value, reason, _dt_to_str(now), _dt_to_str(expires)),
            )
            con.commit()

    def has_negative(self, provider: str, query_type: str, query_value: str) -> bool:
        qhash = hashlib.sha256(query_value.encode("utf-8")).hexdigest()
        with self._connect() as con:
            cur = con.execute(
                "SELECT expires_at FROM negative_cache WHERE provider=? AND query_type=? AND query_hash=?",
                (provider, query_type, qhash),
            )
            row = cur.fetchone()
            if not row:
                return False
            expires_at = datetime.fromisoformat(row[0]) if row[0] else None
            active = not (expires_at and expires_at < _now())
            if active:
                self.stats["negative_cache_hits"] += 1
            return active
