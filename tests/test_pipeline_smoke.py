from pathlib import Path
from gas_screening_mvp.pipeline import run_pipeline


def test_pipeline_smoke(tmp_path):
    out = tmp_path / "out.xlsx"
    sheets = run_pipeline(
        config_path=Path("examples/demo_config.yml"),
        input_csv=Path("examples/sample_input.csv"),
        output_xlsx=out,
    )
    assert out.exists()
    assert "Summary" in sheets
    assert len(sheets["Summary"]) >= 3
    cf4 = [r for r in sheets["Summary"] if r["preferred_name"] == "carbon tetrafluoride"]
    assert cf4
    assert cf4[0]["pfas_flag"] == "yes"


def test_pipeline_can_export_csv_bundle(tmp_path):
    out_dir = tmp_path / "csv_output"
    sheets = run_pipeline(
        config_path=Path("examples/demo_config.yml"),
        input_csv=Path("examples/sample_input.csv"),
        output_xlsx=out_dir,
    )

    assert (out_dir / "summary.csv").exists()
    assert (out_dir / "run_stats.csv").exists()
    assert "Summary" in sheets
