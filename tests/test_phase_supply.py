from gas_screening_mvp.derivation.phase import phase_at_25C
from gas_screening_mvp.derivation.supply_class import derive_supply_class
from gas_screening_mvp.derivation.summary_builder import build_summary_row
from gas_screening_mvp.domain.models import ChemicalIdentity, SelectedProperty


def test_phase_liquid():
    assert phase_at_25C(200.0, 350.0, 500.0) == "liquid"


def test_phase_above_tc():
    assert phase_at_25C(50.0, 100.0, 200.0) == "gas_or_supercritical"


def test_supply_gas():
    assert derive_supply_class("gas", None, None) == "compressed_or_liquefied_gas"


def test_supply_heated():
    assert derive_supply_class("liquid", 0.5, 2.0) == "heated_source_required"


def test_summary_marks_vapor_pressure_above_tc_not_applicable():
    chemical = ChemicalIdentity(
        chemical_id="chem-1",
        candidate_id="cand-1",
        preferred_name="carbon tetrafluoride",
        cas="75-73-0",
        pubchem_cid=None,
        formula="CF4",
        molecular_weight=88.003,
        canonical_smiles="FC(F)(F)F",
        isomeric_smiles="FC(F)(F)F",
        inchi=None,
        inchikey=None,
        identity_status="manual_review_required",
        confidence=0.35,
        source="test",
    )
    selected = [
        SelectedProperty("chem-1", "normal_melting_point", 89.56, None, "K", "selected", "C", "CuratedCsv", "test"),
        SelectedProperty("chem-1", "normal_boiling_point", 145.1, None, "K", "selected", "C", "CuratedCsv", "test"),
        SelectedProperty("chem-1", "critical_temperature", 227.51, None, "K", "selected", "C", "CuratedCsv", "test"),
        SelectedProperty("chem-1", "critical_pressure", 3740000, None, "Pa", "selected", "C", "CuratedCsv", "test"),
        SelectedProperty("chem-1", "vapor_pressure_298.15K", None, None, "Pa", "missing", "Missing", None, "test"),
    ]

    row = build_summary_row(chemical, selected, [], [])

    assert row["phase_25C_1atm"] == "gas_or_supercritical"
    assert row["pvap_25C_status"] == "not_applicable_above_Tc"
