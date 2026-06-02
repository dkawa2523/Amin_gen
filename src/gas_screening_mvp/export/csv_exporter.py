from __future__ import annotations

import csv
import re
from pathlib import Path


def export_csv(sheets: dict[str, list[dict]], output_path: str | Path) -> None:
    output = Path(output_path)
    if output.suffix.lower() == ".csv":
        output.parent.mkdir(parents=True, exist_ok=True)
        _write_csv(output, sheets.get("Summary", []))
        return

    output.mkdir(parents=True, exist_ok=True)
    for sheet_name, rows in sheets.items():
        _write_csv(output / f"{_slugify(sheet_name)}.csv", rows)


def _write_csv(path: Path, rows: list[dict]) -> None:
    headers = _headers(rows)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        if rows:
            writer.writerows(rows)
        else:
            writer.writerow({"message": "no rows"})


def _headers(rows: list[dict]) -> list[str]:
    if not rows:
        return ["message"]
    headers: list[str] = []
    for row in rows:
        for key in row:
            if key not in headers:
                headers.append(key)
    return headers


def _slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower()).strip("_")
    return slug or "sheet"
