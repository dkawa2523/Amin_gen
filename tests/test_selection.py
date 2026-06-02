from gas_screening_mvp.domain.models import PropertyCandidate
from gas_screening_mvp.selection.property_selector import PropertySelector
from gas_screening_mvp.selection.unit_normalizer import normalize_property


def test_degree_celsius_symbol_converts_to_kelvin():
    prop = PropertyCandidate("chem-1", "normal_boiling_point", value_num=25.0, unit="\u00b0C")

    normalized = normalize_property(prop, "K")

    assert normalized.value_num == 298.15
    assert normalized.unit == "K"


def test_curated_csv_priority_wins_by_default():
    selector = PropertySelector()
    selected = selector.select(
        [
            PropertyCandidate("chem-1", "critical_temperature", value_num=508.0, unit="K", source="CoolProp", quality_hint="A"),
            PropertyCandidate("chem-1", "critical_temperature", value_num=510.0, unit="K", source="CuratedCsv", quality_hint="A"),
        ],
        ["critical_temperature"],
    )

    assert selected[0].value_num == 510.0
    assert selected[0].selected_source == "CuratedCsv"


def test_vapor_pressure_selection_marks_outside_valid_range():
    selector = PropertySelector()
    selected = selector.select_by_temperature(
        [
            PropertyCandidate(
                "chem-1",
                "vapor_pressure",
                value_num=1000.0,
                unit="Pa",
                temperature_K=298.15,
                source="thermo.VaporPressure",
                quality_hint="B",
                valid_temperature_min_K=310.0,
                valid_temperature_max_K=360.0,
            )
        ],
        "vapor_pressure",
        298.15,
    )

    assert selected.status == "outside_range"
    assert selected.quality_rank == "D"
