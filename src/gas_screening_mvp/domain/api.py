from __future__ import annotations

from dataclasses import dataclass
from typing import Literal
import hashlib
import json


@dataclass(frozen=True)
class ApiRequest:
    provider: str
    endpoint: str
    method: Literal["GET", "POST"]
    params: dict
    body: dict | None = None
    headers: dict | None = None

    def signature(self) -> str:
        payload = {
            "provider": self.provider,
            "endpoint": self.endpoint,
            "method": self.method,
            "params": self.params,
            "body": self.body,
        }
        text = json.dumps(payload, sort_keys=True, ensure_ascii=False, default=str)
        return hashlib.sha256(text.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class ApiResponse:
    request: ApiRequest
    status_code: int
    text: str
    error_type: str | None = None
