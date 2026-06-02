from __future__ import annotations

import time
from collections import defaultdict


class SyncRateLimiter:
    def __init__(self, rates: dict[str, float]):
        self.rates = rates
        self.last_call = defaultdict(lambda: 0.0)

    def wait(self, provider: str) -> None:
        rps = float(self.rates.get(provider, 1.0))
        if rps <= 0:
            rps = 1.0
        min_interval = 1.0 / rps
        elapsed = time.time() - self.last_call[provider]
        if elapsed < min_interval:
            time.sleep(min_interval - elapsed)
        self.last_call[provider] = time.time()
