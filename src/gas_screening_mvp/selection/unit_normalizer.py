from __future__ import annotations

from dataclasses import replace

from gas_screening_mvp.domain.models import PropertyCandidate


def normalize_property(candidate: PropertyCandidate, target_unit: str | None) -> PropertyCandidate:
    if target_unit is None:
        return candidate
    value = candidate.value_num
    if value is None:
        return candidate

    src = (candidate.unit or "").strip()
    dst = target_unit.strip()
    src_key = _unit_key(src)
    dst_key = _unit_key(dst)
    if src_key == dst_key:
        return replace(candidate, unit=dst)

    new_value = value
    if src_key == "C" and dst_key == "K":
        new_value = value + 273.15
    elif src_key == "K" and dst_key == "C":
        new_value = value - 273.15
    elif src_key == "PA" and dst_key == "KPA":
        new_value = value / 1000.0
    elif src_key == "KPA" and dst_key == "PA":
        new_value = value * 1000.0
    elif src_key == "PA" and dst_key == "MPA":
        new_value = value / 1_000_000.0
    elif src_key == "MPA" and dst_key == "PA":
        new_value = value * 1_000_000.0
    elif src_key == "BAR" and dst_key == "PA":
        new_value = value * 100_000.0
    elif src_key == "PA" and dst_key == "BAR":
        new_value = value / 100_000.0
    else:
        return candidate

    return replace(candidate, value_num=float(new_value), unit=dst)


def _unit_key(unit: str) -> str:
    compact = unit.strip().replace(" ", "")
    if compact in {"C", "°C", "℃", "degC", "DEGC", "celsius", "Celsius"}:
        return "C"
    if compact in {"K", "kelvin", "Kelvin"}:
        return "K"
    pressure = {
        "Pa": "PA",
        "pascal": "PA",
        "Pascal": "PA",
        "kPa": "KPA",
        "KPa": "KPA",
        "kilopascal": "KPA",
        "MPa": "MPA",
        "megapascal": "MPA",
        "bar": "BAR",
    }
    return pressure.get(compact, compact.upper())
