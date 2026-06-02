from __future__ import annotations

import argparse
from pathlib import Path

from gas_screening_mvp.pipeline import run_pipeline


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Semiconductor gas candidate screening MVP")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Run the screening pipeline")
    run.add_argument("--config", default=None, help="Path to YAML config")
    run.add_argument("--input", default=None, help="Optional input CSV with name,smiles,cas,family")
    run.add_argument("--output", required=True, help="Output XLSX path, CSV path, or CSV output directory")
    run.add_argument("--mode", choices=["exploration", "enrichment", "refresh"], default=None, help="Run mode override")
    remote = run.add_mutually_exclusive_group()
    remote.add_argument("--remote", dest="remote", action="store_true", help="Enable configured remote PubChem/PUG-View enrichment")
    remote.add_argument("--no-remote", dest="remote", action="store_false", help="Disable remote PubChem/PUG-View enrichment")
    run.set_defaults(remote=None)
    run.add_argument("--max-candidates", type=int, default=None, help="Override generation.max_total_candidates")
    run.add_argument("--dry-run", action="store_true", help="Plan generation/prefilter/API work without calling remote APIs")

    args = parser.parse_args(argv)
    if args.cmd == "run":
        run_pipeline(
            args.config,
            args.input,
            args.output,
            mode=args.mode,
            remote=args.remote,
            max_candidates=args.max_candidates,
            dry_run=args.dry_run,
        )
        suffix = " dry-run plan" if args.dry_run else ""
        print(f"Wrote{suffix} {Path(args.output).resolve()}")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
