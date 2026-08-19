from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from gmsl_budget.icgem import discover_gfc_downloads, download_gfc


SERIES_PAGES = (
    "https://icgem.gfz.de/sp/01_GRACE/CSR/CSR%20Release%2006",
    "https://icgem.gfz.de/sp/01_GRACE/CSR/CSR%20Release%2006.3%20%28GFO%29",
)
SOURCE_METADATA = {
    "solution_code": "BA01 (unconstrained degree/order 60; BB01 degree/order 96 excluded)",
    "duplicate_month_rule": "retain the BA01 arc whose epoch midpoint is closest to day 15; prefer longer arc on ties",
    "GRACE": {"release": "CSR RL06", "license": "CC BY 4.0"},
    "GRACE-FO": {
        "release": "CSR RL06.3",
        "doi": "10.5067/GFL20-MC063",
        "license": "CC BY 4.0",
    },
}


def main() -> int:
    parser = argparse.ArgumentParser(description="Download official ICGEM CSR DDK1 monthly GFC files.")
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--start", default="2013-11")
    parser.add_argument("--end", default="2024-10")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    specs = discover_gfc_downloads(SERIES_PAGES, args.start, args.end)
    if not specs:
        raise RuntimeError("ICGEM discovery returned no CSR DDK1 files")
    print(f"Discovered {len(specs)} CSR DDK1 epochs from {args.start} through {args.end}")
    for spec in specs:
        print(spec.filename)
    if args.dry_run:
        return 0
    results = [download_gfc(spec, args.destination) for spec in specs]
    manifest = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "start_month": args.start,
        "end_month": args.end,
        "series_pages": list(SERIES_PAGES),
        "source_metadata": SOURCE_METADATA,
        "files": [
            {
                "filename": result.path.name,
                "source_url": result.source_url,
                "source_page": result.source_page,
                "sha256": result.sha256,
                "size_bytes": result.size_bytes,
            }
            for result in results
        ],
    }
    manifest_path = args.destination / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
