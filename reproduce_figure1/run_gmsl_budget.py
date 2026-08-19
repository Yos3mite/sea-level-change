from __future__ import annotations

import argparse

from gmsl_budget.config import PipelineConfig
from gmsl_budget.pipeline import run_pipeline


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the auditable full-ocean GMSL budget workflow")
    parser.add_argument("--config", required=True, help="JSON configuration path")
    args = parser.parse_args()
    output = run_pipeline(PipelineConfig.load(args.config))
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
