from pathlib import Path

from gas_screening_mvp.config_loader import load_config
from gas_screening_mvp.cli import main
from gas_screening_mvp.domain.models import ChemicalIdentity, MoleculeCandidate, NormalizedMolecule, PropertyCandidate, SelectedProperty
from gas_screening_mvp.generation.amine_template import AmineTemplateGenerator
from gas_screening_mvp.generation.local_mutation import LocalMutationGenerator
from gas_screening_mvp.generation.semiconductor_amine import SemiconductorAmineGenerator
from gas_screening_mvp.pipeline import ScreeningPipeline, run_pipeline
from gas_screening_mvp.planning.fetch_planner import FetchPlanner
from gas_screening_mvp.providers.amine_estimates import AmineLocalEstimateProvider
from gas_screening_mvp.selection.property_selector import PropertySelector
from gas_screening_mvp.storage.cache import SqliteApiCache


def test_amine_template_generator_has_no_duplicate_candidates():
    gen = AmineTemplateGenerator(
        max_candidates=60,
        max_carbons=5,
        max_heavy_atoms=12,
        allowed_substituents=["small_alkyl", "fluorinated_alkyl"],
    )

    candidates = list(gen.generate())
    ids = [candidate.candidate_id for candidate in candidates]
    smiles = [candidate.smiles for candidate in candidates]

    assert candidates
    assert len(ids) == len(set(ids))
    assert len(smiles) == len(set(smiles))


def test_amine_template_generator_adds_classification_metadata():
    gen = AmineTemplateGenerator(
        max_candidates=20,
        max_carbons=4,
        max_heavy_atoms=12,
        allowed_substituents=["methyl", "trifluoromethyl"],
    )

    candidates = list(gen.generate())
    fluorinated = [candidate for candidate in candidates if candidate.metadata.get("contains_fluorinated_substituent") == "true"]

    assert {candidate.metadata.get("amine_class") for candidate in candidates} >= {"primary", "secondary"}
    assert fluorinated
    assert all(candidate.metadata.get("amine_substituents") is not None for candidate in candidates if candidate.generation_rule != "cyclic_amine")


def test_semiconductor_amine_generator_adds_process_relevance_metadata():
    candidates = list(SemiconductorAmineGenerator().generate())

    assert candidates
    assert any(candidate.metadata.get("precursor_family") == "aminosilane_precursor" for candidate in candidates)
    assert any("ald_cvd_precursor_candidate" in candidate.metadata.get("semiconductor_process_roles", "") for candidate in candidates)
    assert len({candidate.smiles for candidate in candidates}) == len(candidates)


def test_amine_local_estimate_provider_returns_temperature_and_vapor_pressure():
    chemical = ChemicalIdentity(
        chemical_id="chem-1",
        candidate_id="cand-1",
        preferred_name="generated_tertiary_amine_methyl_methyl_ethyl",
        cas=None,
        pubchem_cid=None,
        formula="C4H11N",
        molecular_weight=73.139,
        canonical_smiles="N(C)(C)CC",
        isomeric_smiles="N(C)(C)CC",
        inchi=None,
        inchikey=None,
        identity_status="manual_review_required",
        confidence=0.35,
        source="test",
    )

    values = AmineLocalEstimateProvider([298.15]).fetch(chemical, ["normal_boiling_point", "vapor_pressure"])

    assert {value.property_name for value in values} == {"normal_boiling_point", "vapor_pressure"}
    assert all(value.source == "LocalAmineEstimate" for value in values)
    assert all(value.is_estimated for value in values)


def test_estimated_property_does_not_conflict_with_curated_value():
    selector = PropertySelector()

    selected = selector.select(
        [
            PropertyCandidate("chem-1", "normal_boiling_point", value_num=362.0, unit="K", source="CuratedCsv", quality_hint="A"),
            PropertyCandidate("chem-1", "normal_boiling_point", value_num=352.0, unit="K", source="LocalAmineEstimate", quality_hint="D", is_estimated=True),
        ],
        ["normal_boiling_point"],
    )

    assert selected[0].value_num == 362.0
    assert selected[0].selected_source == "CuratedCsv"
    assert selected[0].status == "selected"


def test_local_mutation_generator_returns_deduped_candidates():
    seeds = [
        MoleculeCandidate("seed-1", "seed", "amine", input_name="methylamine", smiles="CN", generation_score=1.0),
        MoleculeCandidate("seed-2", "seed", "amine", input_name="ethylamine", smiles="CCN", generation_score=1.0),
    ]
    gen = LocalMutationGenerator(seeds, enabled=True, max_candidates=20)

    candidates = list(gen.generate())
    smiles = [candidate.smiles for candidate in candidates]

    assert candidates
    assert len(smiles) == len(set(smiles))
    assert all(candidate.source == "local_mutation" for candidate in candidates)


def test_formula_only_candidates_do_not_plan_remote_when_remote_disabled(tmp_path):
    cfg = load_config(Path("examples/exploration_config.yml"))
    cfg["providers"]["cache_db"] = str(tmp_path / "cache.sqlite")
    cfg["providers"]["pubchem_enabled"] = False
    cfg["providers"]["pugview_enabled"] = False
    cfg["generation"]["enabled_generators"] = ["FluorocarbonFormulaGenerator"]
    cfg["generation"]["include_formula_only"] = True
    cfg["generation"]["max_total_candidates"] = 25
    pipe = ScreeningPipeline(cfg)

    sheets = pipe.run(dry_run=True)
    stats = {row["metric"]: row["value"] for row in sheets["Run Stats"]}

    assert stats["generated_candidates"] > 10
    assert stats["pubchem_requests"] == 0
    assert stats["planned_api_requests"] == 0
    assert stats["planned_api_candidates"] == 0
    assert sheets["Planned API"] == []


def test_fetch_planner_skips_below_threshold(tmp_path):
    cache = SqliteApiCache(tmp_path / "cache.sqlite")
    planner = FetchPlanner(cache, {"min_identity_api_score": 0.95})
    mol = _normalized("candidate-1", canonical_smiles="FC(F)(F)F")
    selected = [_missing_property("candidate-1")]

    assert planner.plan_for_molecule(mol, selected, remote_enabled=True) == []


def test_fetch_planner_skips_negative_cache(tmp_path):
    cache = SqliteApiCache(tmp_path / "cache.sqlite")
    planner = FetchPlanner(cache, {"min_identity_api_score": 0.1})
    mol = _normalized("candidate-1", canonical_smiles="FC(F)(F)F")
    selected = [_missing_property("candidate-1")]
    cache.put_negative("PubChemPugRest", "smiles", "FC(F)(F)F", "test")

    decisions = planner.plan_for_molecule(mol, selected, remote_enabled=True)

    assert not [decision for decision in decisions if decision.provider == "PubChemPugRest"]


def test_fetch_planner_respects_individual_provider_flags(tmp_path):
    cache = SqliteApiCache(tmp_path / "cache.sqlite")
    planner = FetchPlanner(cache, {"min_identity_api_score": 0.1})
    mol = _normalized("candidate-1", canonical_smiles="FC(F)(F)F")
    selected = [_missing_property("candidate-1")]

    decisions = planner.plan_for_molecule(
        mol,
        selected,
        remote_enabled=True,
        pubchem_enabled=False,
        pugview_enabled=True,
    )

    assert not [decision for decision in decisions if decision.provider == "PubChemPugRest"]


def test_dry_run_completes_without_api(tmp_path):
    out = tmp_path / "dry_run.xlsx"

    sheets = run_pipeline(
        config_path=Path("examples/demo_config.yml"),
        input_csv=Path("examples/sample_input.csv"),
        output_xlsx=out,
        dry_run=True,
    )
    stats = {row["metric"]: row["value"] for row in sheets["Run Stats"]}

    assert out.exists()
    assert stats["pubchem_requests"] == 0
    assert stats["pugview_requests"] == 0
    assert stats["planned_api_requests"] == 0
    assert stats["planned_api_candidates"] == 0
    assert sheets["Summary"] == []


def test_exploration_config_pipeline_runs_without_remote_api(tmp_path):
    cfg_path = Path("examples/exploration_config.yml")

    sheets = run_pipeline(
        config_path=cfg_path,
        input_csv=Path("examples/sample_input.csv"),
        output_xlsx=tmp_path / "exploration.xlsx",
        dry_run=True,
    )
    stats = {row["metric"]: row["value"] for row in sheets["Run Stats"]}

    assert stats["generated_candidates"] > 5
    assert stats["pubchem_requests"] == 0
    assert stats["planned_api_requests"] == 0


def test_pipeline_summary_includes_candidate_annotations(tmp_path):
    sheets = run_pipeline(
        config_path=Path("examples/evaluation_exploration_config.yml"),
        input_csv=Path("examples/sample_input.csv"),
        output_xlsx=tmp_path / "annotated.csv",
        mode="exploration",
        remote=False,
        max_candidates=40,
    )

    amines = [row for row in sheets["Summary"] if row.get("candidate_family") == "amine"]

    assert amines
    assert any(row.get("amine_class") == "tertiary" for row in amines)
    assert any(row.get("amine_class_label") == "tertiary_amine" for row in amines)
    assert all(row.get("reactive_groups") for row in amines)
    assert any(row.get("precursor_family") == "aminosilane_precursor" for row in amines)
    assert any(row.get("semiconductor_process_roles") for row in amines)
    assert any(row.get("tb_C") not in (None, "") for row in amines)
    assert "Candidate Breakdown" in sheets
    assert "Amine Summary" in sheets
    assert "Amine Class Summary" in sheets
    assert any(row["breakdown"] == "amine_class" for row in sheets["Candidate Breakdown"])
    assert any(row["amine_class"] == "tertiary" for row in sheets["Amine Class Summary"])


def test_cli_dry_run_completes_without_api(tmp_path):
    out = tmp_path / "cli_dry_run.xlsx"

    code = main([
        "run",
        "--config", "examples/exploration_config.yml",
        "--input", "examples/sample_input.csv",
        "--output", str(out),
        "--mode", "exploration",
        "--dry-run",
        "--max-candidates", "30",
    ])

    assert code == 0
    assert out.exists()


def _normalized(candidate_id: str, canonical_smiles: str | None = None, formula: str | None = "CF4") -> NormalizedMolecule:
    return NormalizedMolecule(
        candidate_id=candidate_id,
        source="seed",
        family="fluorocarbon",
        input_name="candidate",
        cas=None,
        canonical_smiles=canonical_smiles,
        isomeric_smiles=canonical_smiles,
        standard_inchi=None,
        standard_inchikey=None,
        formula=formula,
        molecular_weight=88.0,
        heavy_atom_count=5,
        element_symbols=("C", "F"),
        structure_status="valid" if canonical_smiles else "manual_review_required",
    )


def _missing_property(chemical_id: str) -> SelectedProperty:
    return SelectedProperty(
        chemical_id=chemical_id,
        property_name="normal_boiling_point",
        value_num=None,
        value_text=None,
        unit="K",
        status="missing",
        quality_rank="Missing",
        selected_source=None,
        selection_reason="test",
    )
