from __future__ import annotations

import math
import re
from collections import Counter
from typing import Iterable

from gas_screening_mvp.domain.models import ChemicalIdentity, PropertyCandidate


class AmineLocalEstimateProvider:
    """Low-confidence local estimates for generated amine candidates.

    These values are meant for offline screening and plotting when no curated
    data are available. They must not be treated as validated design data.
    """

    name = "LocalAmineEstimate"

    def __init__(self, target_temperatures_K: list[float] | None = None):
        self.target_temperatures_K = target_temperatures_K or [298.15, 313.15, 333.15]

    def supports(self, property_name: str) -> bool:
        return property_name in {
            "normal_melting_point",
            "normal_boiling_point",
            "critical_temperature",
            "critical_pressure",
            "vapor_pressure",
        }

    def fetch(self, chemical: ChemicalIdentity, property_names: Iterable[str]) -> list[PropertyCandidate]:
        wanted = set(property_names)
        counts = _formula_counts(chemical.formula or "")
        if not _looks_like_amine_candidate(chemical, counts):
            return []

        mw = chemical.molecular_weight or _estimate_mw(counts)
        if mw is None:
            return []

        boiling_K = _estimate_boiling_point_K(chemical, counts, mw)
        melting_K = _estimate_melting_point_K(chemical, counts, mw, boiling_K)
        critical_K = max(boiling_K + 90.0, _estimate_critical_temperature_K(boiling_K, mw))
        critical_pressure_Pa = _estimate_critical_pressure_Pa(mw, counts)

        out: list[PropertyCandidate] = []
        if "normal_melting_point" in wanted:
            out.append(_candidate(chemical, "normal_melting_point", melting_K, "K"))
        if "normal_boiling_point" in wanted:
            out.append(_candidate(chemical, "normal_boiling_point", boiling_K, "K"))
        if "critical_temperature" in wanted:
            out.append(_candidate(chemical, "critical_temperature", critical_K, "K"))
        if "critical_pressure" in wanted:
            out.append(_candidate(chemical, "critical_pressure", critical_pressure_Pa, "Pa"))
        if "vapor_pressure" in wanted:
            for T in self.target_temperatures_K:
                if T >= critical_K:
                    continue
                p = _estimate_vapor_pressure_Pa(T, boiling_K, mw, counts)
                out.append(_candidate(chemical, "vapor_pressure", p, "Pa", temperature_K=T))
        return out


def _candidate(
    chemical: ChemicalIdentity,
    property_name: str,
    value: float,
    unit: str,
    temperature_K: float | None = None,
) -> PropertyCandidate:
    return PropertyCandidate(
        chemical_id=chemical.chemical_id,
        property_name=property_name,
        value_num=round(float(value), 6),
        unit=unit,
        temperature_K=temperature_K,
        source=AmineLocalEstimateProvider.name,
        method="rough_group_heuristic",
        is_estimated=True,
        quality_hint="D",
        reference="local screening estimate; replace with curated or validated property data before design use",
    )


def _looks_like_amine_candidate(chemical: ChemicalIdentity, counts: Counter[str]) -> bool:
    name = (chemical.preferred_name or "").lower()
    smiles = chemical.canonical_smiles or chemical.isomeric_smiles or ""
    if counts.get("N", 0) <= 0:
        return False
    if name in {"nitrogen", "nitrogen trifluoride", "nitrous oxide"}:
        return False
    if name == "ammonia":
        return True
    if any(token in name for token in ("amine", "amino", "amide", "azane", "piper", "pyrrol", "morpholine")):
        return True
    return bool(counts.get("C", 0) or counts.get("Si", 0) or counts.get("B", 0) or counts.get("Ti", 0) or "N(" in smiles)


def _estimate_boiling_point_K(chemical: ChemicalIdentity, counts: Counter[str], mw: float) -> float:
    name = (chemical.preferred_name or "").lower()
    if name == "ammonia":
        return 239.8
    c = counts.get("C", 0)
    n = counts.get("N", 0)
    f = counts.get("F", 0)
    si = counts.get("Si", 0)
    b = counts.get("B", 0)
    metal = counts.get("Ti", 0) + counts.get("W", 0) + counts.get("Mo", 0)
    ring_bonus = 12.0 if any(token in name for token in ("piper", "pyrrol", "morpholine", "aziridine", "azetidine")) else 0.0
    tb = 210.0 + 1.30 * mw + 8.0 * max(0, n - 1) + 4.0 * f + 16.0 * si + 10.0 * b + 28.0 * metal + ring_bonus
    if c <= 1 and n == 1 and not f and not si:
        tb -= 12.0
    return min(max(tb, 210.0), 620.0)


def _estimate_melting_point_K(chemical: ChemicalIdentity, counts: Counter[str], mw: float, boiling_K: float) -> float:
    name = (chemical.preferred_name or "").lower()
    if name == "ammonia":
        return 195.4
    symmetry_bonus = 18.0 if counts.get("Si", 0) or counts.get("B", 0) or counts.get("Ti", 0) else 0.0
    fluorine_bonus = min(counts.get("F", 0) * 1.5, 18.0)
    tm = boiling_K - 165.0 + 0.08 * mw + symmetry_bonus + fluorine_bonus
    return min(max(tm, 95.0), boiling_K - 10.0)


def _estimate_critical_temperature_K(boiling_K: float, mw: float) -> float:
    return boiling_K + 145.0 + 0.22 * mw


def _estimate_critical_pressure_Pa(mw: float, counts: Counter[str]) -> float:
    hetero_bonus = 120000.0 * (counts.get("N", 0) + counts.get("O", 0))
    heavy_penalty = 8500.0 * mw
    return min(max(5_200_000.0 - heavy_penalty + hetero_bonus, 1_200_000.0), 6_800_000.0)


def _estimate_vapor_pressure_Pa(T_K: float, boiling_K: float, mw: float, counts: Counter[str]) -> float:
    delta_h_vap_J_mol = (30.0 + 0.09 * mw + 1.2 * counts.get("N", 0) + 0.4 * counts.get("F", 0)) * 1000.0
    ln_p = math.log(101325.0) - delta_h_vap_J_mol / 8.314462618 * (1.0 / T_K - 1.0 / boiling_K)
    return min(max(math.exp(ln_p), 1e-6), 8_000_000.0)


def _formula_counts(formula: str) -> Counter[str]:
    counts: Counter[str] = Counter()
    for symbol, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[symbol] += int(count or "1")
    return counts


def _estimate_mw(counts: Counter[str]) -> float | None:
    weights = {
        "H": 1.008,
        "B": 10.81,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "F": 18.998,
        "Si": 28.085,
        "Ti": 47.867,
    }
    total = 0.0
    for symbol, count in counts.items():
        weight = weights.get(symbol)
        if weight is None:
            return None
        total += weight * count
    return total
