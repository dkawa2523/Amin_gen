from pathlib import Path

from gas_screening_mvp.config_loader import load_config
from gas_screening_mvp.domain.models import NormalizedMolecule
from gas_screening_mvp.generation.registry import generate_candidates
from gas_screening_mvp.pipeline import ScreeningPipeline
from gas_screening_mvp.providers.pubchem_pugrest import PubChemPugRestProvider
from gas_screening_mvp.storage.cache import SqliteApiCache


class DummyResponse:
    def __init__(self, status_code: int, text: str):
        self.status_code = status_code
        self.text = text


def test_pubchem_expansion_generator_uses_mocked_api_and_cache(tmp_path, monkeypatch):
    calls = []

    def fake_post(endpoint, params=None, data=None, timeout=None):
        calls.append((endpoint, params, data, timeout))
        return DummyResponse(200, "123\n456\n789\n")

    import gas_screening_mvp.generation.pubchem_expansion as expansion_module

    monkeypatch.setattr(expansion_module.requests, "post", fake_post)
    cache = SqliteApiCache(tmp_path / "cache.sqlite")
    cfg = {
        "max_total_candidates": 10,
        "enabled_generators": ["SeedListGenerator", "PubChemExpansionGenerator"],
        "pubchem_expansion_enabled": True,
        "pubchem_similarity_threshold": 91,
        "max_cids_per_seed": 2,
        "pubchem_expansion_rps": 1000.0,
    }

    first = generate_candidates(
        cfg,
        input_csv=Path("examples/sample_input.csv"),
        cache=cache,
        mode="exploration",
        remote_enabled=True,
        dry_run=False,
    )
    second = generate_candidates(
        cfg,
        input_csv=Path("examples/sample_input.csv"),
        cache=cache,
        mode="exploration",
        remote_enabled=True,
        dry_run=False,
    )

    expanded = [candidate for candidate in first if candidate.source == "pubchem_similarity"]
    assert expanded
    assert all(candidate.metadata.get("pubchem_cid") for candidate in expanded)
    assert len(calls) == 5
    assert [candidate.candidate_id for candidate in first] == [candidate.candidate_id for candidate in second]


def test_pubchem_expansion_not_connected_outside_exploration(tmp_path, monkeypatch):
    def fail_post(*args, **kwargs):
        raise AssertionError("remote API should not be called")

    import gas_screening_mvp.generation.pubchem_expansion as expansion_module

    monkeypatch.setattr(expansion_module.requests, "post", fail_post)
    cfg = {
        "max_total_candidates": 10,
        "enabled_generators": ["SeedListGenerator", "PubChemExpansionGenerator"],
        "pubchem_expansion_enabled": True,
        "max_cids_per_seed": 2,
    }

    candidates = generate_candidates(
        cfg,
        input_csv=Path("examples/sample_input.csv"),
        cache=SqliteApiCache(tmp_path / "cache.sqlite"),
        mode="enrichment",
        remote_enabled=True,
        dry_run=False,
    )

    assert candidates
    assert not [candidate for candidate in candidates if candidate.source == "pubchem_similarity"]


def test_formula_search_resolves_formula_only_candidate_with_cache(tmp_path, monkeypatch):
    calls = []

    def fake_get(endpoint, params=None, timeout=None):
        calls.append((endpoint, params, timeout))
        if "fastformula" in endpoint:
            return DummyResponse(200, "123\n999\n")
        if "/property/" in endpoint:
            return DummyResponse(
                200,
                "CID,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName\n"
                "123,CF4,88.003,FC(F)(F)F,FC(F)(F)F,InChI=1S/CF4,TXEYQDLBPFQVAA-UHFFFAOYSA-N,tetrafluoromethane\n",
            )
        raise AssertionError(f"unexpected endpoint: {endpoint}")

    import gas_screening_mvp.providers.pubchem_pugrest as pugrest_module

    monkeypatch.setattr(pugrest_module.requests, "get", fake_get)
    provider = PubChemPugRestProvider(SqliteApiCache(tmp_path / "cache.sqlite"), enabled=True, rps=1000.0)
    molecule = _formula_only_molecule("CF4")

    first = provider.resolve_formula(
        molecule,
        max_cids_per_formula=2,
        allowed_elements={"C", "F"},
        max_heavy_atoms=5,
        max_molecular_weight=100.0,
    )
    second = provider.resolve_formula(
        molecule,
        max_cids_per_formula=2,
        allowed_elements={"C", "F"},
        max_heavy_atoms=5,
        max_molecular_weight=100.0,
    )

    assert first
    assert first[0].canonical_smiles == "FC(F)(F)F"
    assert first[0].inchikey == "TXEYQDLBPFQVAA-UHFFFAOYSA-N"
    assert [hit.pubchem_cid for hit in first] == [hit.pubchem_cid for hit in second]
    assert len(calls) == 2


def test_pipeline_formula_search_connection_uses_mocked_api(tmp_path, monkeypatch):
    input_csv = tmp_path / "formula_only.csv"
    input_csv.write_text("name,smiles,formula,cas,family\ncf4_formula,,CF4,,fluorocarbon\n", encoding="utf-8")
    calls = []

    def fake_get(endpoint, params=None, timeout=None):
        calls.append((endpoint, params, timeout))
        if "fastformula" in endpoint:
            return DummyResponse(200, "123\n")
        if "/property/" in endpoint:
            return DummyResponse(
                200,
                "CID,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName\n"
                "123,CF4,88.003,FC(F)(F)F,FC(F)(F)F,InChI=1S/CF4,TXEYQDLBPFQVAA-UHFFFAOYSA-N,tetrafluoromethane\n",
            )
        raise AssertionError(f"unexpected endpoint: {endpoint}")

    import gas_screening_mvp.providers.pubchem_pugrest as pugrest_module

    monkeypatch.setattr(pugrest_module.requests, "get", fake_get)
    cfg = load_config()
    cfg["mode"] = "exploration"
    cfg["generation"]["enabled_generators"] = ["SeedListGenerator"]
    cfg["generation"]["formula_search_enabled"] = True
    cfg["generation"]["max_cids_per_formula"] = 1
    cfg["prefilter"]["allow_formula_only"] = True
    cfg["providers"]["cache_db"] = str(tmp_path / "pipeline_cache.sqlite")
    cfg["providers"]["pubchem_enabled"] = True
    cfg["providers"]["pugview_enabled"] = False
    cfg["providers"]["max_api_candidates_per_run"] = 1
    cfg["fetch_planning"]["min_identity_api_score"] = 0.1

    sheets = ScreeningPipeline(cfg).run(input_csv=input_csv, dry_run=False)

    assert calls
    assert sheets["Summary"]
    assert sheets["Summary"][0]["pubchem_cid"] == 123
    stats = {row["metric"]: row["value"] for row in sheets["Run Stats"]}
    assert stats["pubchem_requests"] == 2


def test_pipeline_uses_local_properties_before_remote_planning(tmp_path, monkeypatch):
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("name,smiles,cas,family\ntriethylamine,CCN(CC)CC,121-44-8,amine\n", encoding="utf-8")
    curated = tmp_path / "curated.csv"
    curated.write_text(
        "key_type,key_value,property_name,value_num,value_text,unit,temperature_K,pressure_Pa,source,method,is_estimated,quality_hint,reference,valid_temperature_min_K,valid_temperature_max_K\n"
        "name,triethylamine,normal_melting_point,158.45,,K,,,CuratedCsv,test,false,A,,,\n"
        "name,triethylamine,normal_boiling_point,362.0,,K,,,CuratedCsv,test,false,A,,,\n"
        "name,triethylamine,critical_temperature,535.0,,K,,,CuratedCsv,test,false,A,,,\n"
        "name,triethylamine,critical_pressure,3040000,,Pa,,,CuratedCsv,test,false,A,,,\n"
        "name,triethylamine,gwp100_ar6,0,,kg_CO2e_per_kg,,,CuratedCsv,test,false,A,,,\n"
        "name,triethylamine,vapor_pressure,7200,,Pa,298.15,,CuratedCsv,test,false,A,,250,400\n"
        "name,triethylamine,vapor_pressure,14000,,Pa,313.15,,CuratedCsv,test,false,A,,250,400\n"
        "name,triethylamine,vapor_pressure,33000,,Pa,333.15,,CuratedCsv,test,false,A,,250,400\n",
        encoding="utf-8",
    )

    def fail_remote(*args, **kwargs):
        raise AssertionError("PubChem should not be called when local data is complete")

    import gas_screening_mvp.providers.pubchem_pugrest as pugrest_module

    monkeypatch.setattr(pugrest_module.requests, "get", fail_remote)
    monkeypatch.setattr(pugrest_module.requests, "post", fail_remote)
    cfg = load_config()
    cfg["providers"]["cache_db"] = str(tmp_path / "cache.sqlite")
    cfg["providers"]["pubchem_enabled"] = True
    cfg["providers"]["pugview_enabled"] = False
    cfg["providers"]["curated_properties_csv"] = str(curated)
    cfg["providers"]["max_api_candidates_per_run"] = 10
    cfg["fetch_planning"]["min_identity_api_score"] = 0.1

    sheets = ScreeningPipeline(cfg).run(input_csv=input_csv)
    stats = {row["metric"]: row["value"] for row in sheets["Run Stats"]}

    assert sheets["Summary"]
    assert sheets["Planned API"] == []
    assert stats["planned_api_requests"] == 0
    assert stats["pubchem_requests"] == 0


def test_pipeline_reports_cache_hits_for_remote_identity(tmp_path, monkeypatch):
    input_csv = tmp_path / "formula_only.csv"
    input_csv.write_text("name,smiles,formula,cas,family\ncf4_formula,,CF4,,fluorocarbon\n", encoding="utf-8")
    calls = []

    def fake_get(endpoint, params=None, timeout=None):
        calls.append((endpoint, params, timeout))
        if "fastformula" in endpoint:
            return DummyResponse(200, "123\n")
        if "/property/" in endpoint:
            return DummyResponse(
                200,
                "CID,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName\n"
                "123,CF4,88.003,FC(F)(F)F,FC(F)(F)F,InChI=1S/CF4,TXEYQDLBPFQVAA-UHFFFAOYSA-N,tetrafluoromethane\n",
            )
        raise AssertionError(f"unexpected endpoint: {endpoint}")

    import gas_screening_mvp.providers.pubchem_pugrest as pugrest_module

    monkeypatch.setattr(pugrest_module.requests, "get", fake_get)
    cfg = load_config()
    cfg["mode"] = "exploration"
    cfg["generation"]["enabled_generators"] = ["SeedListGenerator"]
    cfg["generation"]["formula_search_enabled"] = True
    cfg["prefilter"]["allow_formula_only"] = True
    cfg["providers"]["cache_db"] = str(tmp_path / "cache.sqlite")
    cfg["providers"]["pubchem_enabled"] = True
    cfg["providers"]["pugview_enabled"] = False
    cfg["providers"]["max_api_candidates_per_run"] = 1
    cfg["fetch_planning"]["min_identity_api_score"] = 0.1

    ScreeningPipeline(cfg).run(input_csv=input_csv)
    sheets = ScreeningPipeline(cfg).run(input_csv=input_csv)
    stats = {row["metric"]: row["value"] for row in sheets["Run Stats"]}

    assert len(calls) == 2
    assert stats["cache_hits"] >= 2
    assert stats["pubchem_requests"] == 0


def test_pipeline_reports_negative_cache_hits(tmp_path):
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("name,smiles,cas,family\ncf4,FC(F)(F)F,,fluorocarbon\n", encoding="utf-8")
    cache = SqliteApiCache(tmp_path / "cache.sqlite")
    cache.put_negative("PubChemPugRest", "smiles", "FC(F)(F)F", "test")
    cfg = load_config()
    cfg["providers"]["cache_db"] = str(tmp_path / "cache.sqlite")
    cfg["providers"]["pubchem_enabled"] = True
    cfg["providers"]["pugview_enabled"] = False
    cfg["providers"]["max_api_candidates_per_run"] = 1
    cfg["fetch_planning"]["min_identity_api_score"] = 0.1

    sheets = ScreeningPipeline(cfg).run(input_csv=input_csv)
    stats = {row["metric"]: row["value"] for row in sheets["Run Stats"]}

    assert stats["negative_cache_hits"] >= 1
    assert not [row for row in sheets["Planned API"] if row["provider"] == "PubChemPugRest"]


def test_formula_search_negative_cache_skips_api(tmp_path, monkeypatch):
    def fail_get(*args, **kwargs):
        raise AssertionError("negative cache should prevent remote API")

    import gas_screening_mvp.providers.pubchem_pugrest as pugrest_module

    monkeypatch.setattr(pugrest_module.requests, "get", fail_get)
    cache = SqliteApiCache(tmp_path / "cache.sqlite")
    cache.put_negative("PubChemPugRest", "formula", "CF4", "test")
    provider = PubChemPugRestProvider(cache, enabled=True)

    assert provider.resolve_formula(_formula_only_molecule("CF4")) == []


def _formula_only_molecule(formula: str) -> NormalizedMolecule:
    return NormalizedMolecule(
        candidate_id="formula-only",
        source="formula_generation",
        family="fluorocarbon",
        input_name=formula,
        cas=None,
        canonical_smiles=None,
        isomeric_smiles=None,
        standard_inchi=None,
        standard_inchikey=None,
        formula=formula,
        molecular_weight=None,
        heavy_atom_count=None,
        element_symbols=(),
        structure_status="manual_review_required",
    )
