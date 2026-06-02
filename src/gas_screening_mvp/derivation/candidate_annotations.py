from __future__ import annotations

from collections import Counter

from gas_screening_mvp.domain.models import ClassificationRecord, NormalizedMolecule


ANNOTATION_SOURCE = "CandidateAnnotation"

AMINE_CLASS_LABELS = {
    "primary": "primary_amine",
    "secondary": "secondary_amine",
    "tertiary": "tertiary_amine",
    "cyclic": "cyclic_amine",
    "diamine": "diamine",
    "polyamine": "polyamine",
    "amino_silane": "amino_silane",
    "silylamine": "silylamine",
    "boron_amide": "boron_amide",
    "metal_amide": "metal_amide",
    "inorganic_nitrogen_source": "inorganic_nitrogen_source",
    "aromatic_heterocycle": "aromatic_nitrogen_heterocycle",
    "unknown": "unknown_amine",
}

AMINE_CLASS_REACTIVE_GROUPS = {
    "primary": "Primary amines / basic nitrogen compounds",
    "secondary": "Secondary amines / basic nitrogen compounds",
    "tertiary": "Tertiary amines / basic nitrogen compounds",
    "cyclic": "Cyclic amines / basic nitrogen compounds",
    "diamine": "Diamines / chelating nitrogen compounds",
    "polyamine": "Polyamines / multidentate nitrogen compounds",
    "amino_silane": "Aminosilanes / silicon precursor amines",
    "silylamine": "Silylamines / surface modifier amines",
    "boron_amide": "Boron amides / boron precursor amines",
    "metal_amide": "Metal amides / ALD-CVD precursor amides",
    "inorganic_nitrogen_source": "Inorganic nitrogen source",
    "aromatic_heterocycle": "Aromatic nitrogen heterocycles",
    "unknown": "Amines / basic nitrogen compounds",
}

SEED_AMINE_ANNOTATIONS = {
    "methylamine": ("primary", ["methyl"]),
    "ethylamine": ("primary", ["ethyl"]),
    "tert-butylamine": ("primary", ["tert-butyl"]),
    "dimethylamine": ("secondary", ["methyl", "methyl"]),
    "diethylamine": ("secondary", ["ethyl", "ethyl"]),
    "trimethylamine": ("tertiary", ["methyl", "methyl", "methyl"]),
    "triethylamine": ("tertiary", ["ethyl", "ethyl", "ethyl"]),
    "morpholine": ("cyclic", []),
    "piperidine": ("cyclic", []),
    "pyridine": ("aromatic_heterocycle", []),
}


def candidate_annotation_row(mol: NormalizedMolecule) -> dict[str, object]:
    metadata = mol.metadata or {}
    amine_class, substituents = _amine_annotation(mol)
    fluorinated_substituent_count = _int_or_zero(metadata.get("fluorinated_substituent_count"))
    if substituents and not fluorinated_substituent_count:
        fluorinated_substituent_count = sum(1 for item in substituents if "fluoro" in item)
    fluorinated_candidate = _contains_element(mol, "F") and mol.family == "amine"
    unsaturated = _bool_or_inferred(metadata.get("contains_unsaturated_substituent"), substituents, {"vinyl", "allyl"})
    cyclic_substituent = _bool_or_inferred(metadata.get("contains_cyclic_substituent"), substituents, {"cyclopropyl"})

    amine_detail = ""
    if amine_class:
        prefix = "fluorinated_" if fluorinated_candidate else ""
        amine_detail = f"{prefix}{amine_class}_amine"

    return {
        "candidate_id": mol.candidate_id,
        "candidate_source": mol.source,
        "candidate_family": mol.family,
        "generation_rule": metadata.get("generation_rule", ""),
        "parent_candidate_id": metadata.get("parent_candidate_id", ""),
        "structure_status": mol.structure_status,
        "structure_status_reason": mol.status_reason,
        "heavy_atom_count": mol.heavy_atom_count,
        "element_symbols": ";".join(mol.element_symbols),
        "amine_class": amine_class,
        "amine_class_label": AMINE_CLASS_LABELS.get(amine_class, "unknown_amine"),
        "amine_detail": amine_detail,
        "amine_substituents": "; ".join(substituents),
        "amine_substituent_profile": _substituent_profile(substituents, amine_class),
        "amine_substituent_count": len(substituents) if amine_class else None,
        "fluorinated_amine": _yes_no(fluorinated_candidate) if mol.family == "amine" else "",
        "fluorinated_substituent_count": fluorinated_substituent_count if mol.family == "amine" else None,
        "amine_fluorination_level": _fluorination_level(fluorinated_candidate, fluorinated_substituent_count) if mol.family == "amine" else "",
        "unsaturated_amine": _yes_no(unsaturated) if mol.family == "amine" else "",
        "cyclic_amine_or_substituent": _yes_no(amine_class == "cyclic" or cyclic_substituent) if mol.family == "amine" else "",
        "ring_name": metadata.get("ring_name", ""),
        "precursor_family": metadata.get("precursor_family", ""),
        "semiconductor_process_roles": metadata.get("semiconductor_process_roles", ""),
        "semiconductor_relevance_basis": metadata.get("semiconductor_relevance_basis", ""),
        "candidate_scope": metadata.get("candidate_scope", ""),
    }


def candidate_classification_records(mol: NormalizedMolecule, chemical_id: str) -> list[ClassificationRecord]:
    annotation = candidate_annotation_row(mol)
    if annotation.get("candidate_family") != "amine":
        return []

    amine_class = str(annotation.get("amine_class") or "unknown")
    out = [
        ClassificationRecord(
            chemical_id,
            "reactive_group",
            AMINE_CLASS_REACTIVE_GROUPS.get(amine_class, AMINE_CLASS_REACTIVE_GROUPS["unknown"]),
            "local generation metadata amine class",
            ANNOTATION_SOURCE,
            confidence="C",
        ),
        ClassificationRecord(
            chemical_id,
            "reactivity_flag",
            "acid_base_reactive",
            "amine/basic nitrogen class",
            ANNOTATION_SOURCE,
            confidence="C",
        ),
    ]
    if "C" in mol.element_symbols:
        out.append(ClassificationRecord(
            chemical_id,
            "reactivity_flag",
            "combustible_or_flammable_review",
            "organic amine screening rule",
            ANNOTATION_SOURCE,
            confidence="D",
        ))
    if annotation.get("fluorinated_amine") == "yes":
        out.append(ClassificationRecord(
            chemical_id,
            "reactivity_flag",
            "fluorinated_amine_pfas_persistence_review",
            "fluorinated amine candidate annotation",
            ANNOTATION_SOURCE,
            confidence="C",
        ))
    if annotation.get("unsaturated_amine") == "yes":
        out.append(ClassificationRecord(
            chemical_id,
            "reactivity_flag",
            "unsaturated_amine_reactivity_review",
            "vinyl/allyl amine candidate annotation",
            ANNOTATION_SOURCE,
            confidence="C",
        ))
    process_roles = annotation.get("semiconductor_process_roles")
    if process_roles:
        out.append(ClassificationRecord(
            chemical_id,
            "reactivity_flag",
            "semiconductor_process_candidate_review",
            str(process_roles),
            ANNOTATION_SOURCE,
            confidence="C",
        ))
    return out


def amine_summary_rows(summary_rows: list[dict]) -> list[dict]:
    columns = [
        "preferred_name",
        "formula",
        "molecular_weight",
        "amine_class_label",
        "amine_class",
        "amine_detail",
        "amine_substituents",
        "amine_substituent_profile",
        "amine_substituent_count",
        "fluorinated_amine",
        "fluorinated_substituent_count",
        "amine_fluorination_level",
        "unsaturated_amine",
        "cyclic_amine_or_substituent",
        "ring_name",
        "precursor_family",
        "semiconductor_process_roles",
        "semiconductor_relevance_basis",
        "gwp100_ar6",
        "gwp100_ar6_status",
        "gwp100_ar6_source",
        "pfas_flag",
        "persistence_screen",
        "reactive_groups",
        "reactivity_flags",
        "phase_25C_1atm",
        "supply_class",
        "tm_C",
        "tb_C",
        "tc_C",
        "pc_MPa",
        "pvap_25C_kPa",
        "pvap_25C_status",
        "pvap_40C_kPa",
        "pvap_40C_status",
        "pvap_60C_kPa",
        "pvap_60C_status",
        "kinetics_coverage",
        "data_quality",
        "review_required",
        "identity_status",
        "candidate_source",
        "generation_rule",
        "structure_status",
        "candidate_id",
    ]
    return [
        {column: row.get(column) for column in columns}
        for row in summary_rows
        if row.get("candidate_family") == "amine"
    ]


def amine_class_summary_rows(summary_rows: list[dict]) -> list[dict]:
    amines = [row for row in summary_rows if row.get("candidate_family") == "amine"]
    classes = sorted(
        {row.get("amine_class") or "unknown" for row in amines},
        key=lambda value: {
            "primary": 1,
            "secondary": 2,
            "tertiary": 3,
            "cyclic": 4,
            "diamine": 5,
            "polyamine": 6,
            "amino_silane": 7,
            "silylamine": 8,
            "boron_amide": 9,
            "metal_amide": 10,
            "inorganic_nitrogen_source": 11,
        }.get(value, 99),
    )
    out = []
    for amine_class in classes:
        rows = [row for row in amines if (row.get("amine_class") or "unknown") == amine_class]
        out.append({
            "amine_class": amine_class,
            "amine_class_label": AMINE_CLASS_LABELS.get(amine_class, "unknown_amine"),
            "count": len(rows),
            "fluorinated_count": _count_value(rows, "fluorinated_amine", "yes"),
            "non_fluorinated_count": _count_value(rows, "fluorinated_amine", "no"),
            "pfas_yes_count": _count_value(rows, "pfas_flag", "yes"),
            "pfas_possible_count": _count_value(rows, "pfas_flag", "possible"),
            "pfas_no_count": _count_value(rows, "pfas_flag", "no"),
            "likely_persistent_count": _count_value(rows, "persistence_screen", "likely_persistent"),
            "reactive_group_assigned_count": sum(1 for row in rows if row.get("reactive_groups")),
            "temperature_property_count": sum(1 for row in rows if all(_present(row.get(field)) for field in ["tm_C", "tb_C", "tc_C", "pc_MPa"])),
            "vapor_pressure_count": sum(1 for row in rows if any(_present(row.get(field)) for field in ["pvap_25C_kPa", "pvap_40C_kPa", "pvap_60C_kPa"])),
            "phase_known_count": sum(1 for row in rows if row.get("phase_25C_1atm") not in {"", None, "unknown"}),
            "review_required_count": _count_value(rows, "review_required", True) + _count_value(rows, "review_required", "True"),
        })
    return out


def candidate_breakdown_rows(summary_rows: list[dict]) -> list[dict]:
    rows: list[dict] = []
    rows.extend(_count_rows(summary_rows, "candidate_family", "candidate_family"))
    rows.extend(_count_rows(summary_rows, "candidate_source", "candidate_source"))
    rows.extend(_count_rows(summary_rows, "generation_rule", "generation_rule"))
    rows.extend(_count_rows(summary_rows, "phase_25C_1atm", "phase_25C_1atm"))
    rows.extend(_count_rows(summary_rows, "supply_class", "supply_class"))
    rows.extend(_count_rows(summary_rows, "data_quality", "data_quality"))

    amines = [row for row in summary_rows if row.get("candidate_family") == "amine"]
    rows.extend(_count_rows(amines, "amine_class", "amine_class"))
    rows.extend(_count_rows(amines, "amine_class_label", "amine_class_label"))
    rows.extend(_count_rows(amines, "amine_substituent_profile", "amine_substituent_profile"))
    rows.extend(_count_rows(amines, "amine_fluorination_level", "amine_fluorination_level"))
    rows.extend(_count_rows(amines, "precursor_family", "amine_precursor_family"))
    rows.extend(_count_rows(amines, "candidate_scope", "candidate_scope"))
    rows.extend(_count_rows(amines, "fluorinated_amine", "fluorinated_amine"))
    rows.extend(_count_rows(amines, "pfas_flag", "amine_pfas_flag"))
    rows.extend(_count_pair_rows(amines, "amine_class", "pfas_flag", "amine_class_by_pfas_flag"))
    rows.extend(_count_pair_rows(amines, "amine_class", "fluorinated_amine", "amine_class_by_fluorination"))
    return rows


def _amine_annotation(mol: NormalizedMolecule) -> tuple[str, list[str]]:
    if mol.family != "amine":
        return "", []
    metadata = mol.metadata or {}
    metadata_class = metadata.get("amine_class")
    metadata_substituents = _split(metadata.get("amine_substituents"))
    if metadata_class:
        return metadata_class, metadata_substituents

    name = (mol.input_name or "").lower()
    if name in SEED_AMINE_ANNOTATIONS:
        return SEED_AMINE_ANNOTATIONS[name]
    if name.startswith("generated_primary_amine_"):
        return "primary", _split_generated_name(name, "generated_primary_amine_")
    if name.startswith("generated_secondary_amine_"):
        return "secondary", _split_generated_name(name, "generated_secondary_amine_")
    if name.startswith("generated_tertiary_amine_"):
        return "tertiary", _split_generated_name(name, "generated_tertiary_amine_")
    if name.startswith("generated_cyclic_amine_"):
        return "cyclic", []
    return "unknown", []


def _split_generated_name(name: str, prefix: str) -> list[str]:
    suffix = name.removeprefix(prefix)
    return [item for item in suffix.split("_") if item]


def _split(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.replace(",", ";").split(";") if item.strip()]


def _substituent_profile(substituents: list[str], amine_class: str) -> str:
    if amine_class == "cyclic":
        return "cyclic_ring"
    if not substituents:
        return ""
    tags = []
    if any("fluoro" in item for item in substituents):
        tags.append("fluorinated_alkyl")
    if any(item in {"vinyl", "allyl"} for item in substituents):
        tags.append("unsaturated")
    if any(item == "cyclopropyl" for item in substituents):
        tags.append("cyclic_alkyl")
    if any("fluoro" not in item and item not in {"vinyl", "allyl", "cyclopropyl"} for item in substituents):
        tags.append("alkyl")
    return "+".join(tags) if tags else "other"


def _fluorination_level(fluorinated_candidate: bool, fluorinated_substituent_count: int) -> str:
    if not fluorinated_candidate:
        return "none"
    if fluorinated_substituent_count <= 1:
        return "single_fluorinated_substituent"
    return "multiple_fluorinated_substituents"


def _contains_element(mol: NormalizedMolecule, element: str) -> bool:
    if element in mol.element_symbols:
        return True
    return element in (mol.formula or "")


def _yes_no(value: bool) -> str:
    return "yes" if value else "no"


def _int_or_zero(value: object) -> int:
    try:
        return int(value)
    except Exception:
        return 0


def _bool_or_inferred(value: str | None, substituents: list[str], markers: set[str]) -> bool:
    if value is not None and value != "":
        return str(value).strip().lower() in {"1", "true", "yes", "y"}
    return any(item in markers for item in substituents)


def _count_value(rows: list[dict], field: str, value: object) -> int:
    return sum(1 for row in rows if row.get(field) == value)


def _present(value: object) -> bool:
    return value not in {None, "", "missing", "unknown", "not_checked"}


def _count_rows(rows: list[dict], field: str, group: str) -> list[dict]:
    counts = Counter((row.get(field) or "unknown") for row in rows)
    return [
        {"breakdown": group, "value": value, "secondary_value": "", "count": count}
        for value, count in sorted(counts.items(), key=lambda item: (-item[1], str(item[0])))
    ]


def _count_pair_rows(rows: list[dict], field_a: str, field_b: str, group: str) -> list[dict]:
    counts = Counter((row.get(field_a) or "unknown", row.get(field_b) or "unknown") for row in rows)
    return [
        {"breakdown": group, "value": value, "secondary_value": secondary, "count": count}
        for (value, secondary), count in sorted(counts.items(), key=lambda item: (-item[1], str(item[0])))
    ]
