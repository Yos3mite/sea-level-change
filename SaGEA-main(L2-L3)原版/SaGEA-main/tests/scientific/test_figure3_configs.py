import json
from pathlib import Path

import numpy as np

from pysrc.reference_products.build_figure3_regional_tws import target_grid_coordinates


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _load(name: str) -> dict:
    return json.loads((PROJECT_ROOT / "config" / name).read_text(encoding="utf-8"))


def test_paper_and_custom_configs_share_the_scientific_regional_contract():
    """Catch custom-L3 comparisons using a different grid, mask, events, or ocean area."""
    paper = _load("figure03_paper_mascon.json")
    custom = _load("figure03_custom_l3.template.json")

    for key in ("target_grid", "mask", "integration", "time", "events", "corrections"):
        assert paper[key] == custom[key]
    assert paper["mode"] == "paper_mascon"
    assert custom["mode"] == "custom_l3"
    assert set(paper["inputs"]["centers"]) == {"CSR", "JPL", "GSFC"}
    assert all(
        "reconstruction" in center
        for center in paper["inputs"]["centers"].values()
    )
    assert custom["inputs"]["custom_l3"]["variables"] == {
        "time": "time",
        "lat": "lat",
        "lon": "lon",
        "field": "field",
        "land_mask": "land_mask",
        "valid_month": "valid_month",
    }
    assert custom["inputs"]["custom_l3"]["metadata"]["gia_corrected"] is False
    assert paper["corrections"] == {"apply_gia": False, "apply_obd": False}


def test_paper_config_contains_only_the_eight_published_numeric_references():
    """Catch plot digitization being mixed with values explicitly reported by Jin et al."""
    paper = _load("figure03_paper_mascon.json")
    references = paper["paper_references"]

    assert references == {
        "el_nino_2014_2016": {
            "total": -6.37,
            "africa": -1.72,
            "north_america": 0.82,
            "south_america": -3.25,
        },
        "el_nino_2023_2024": {
            "total": -4.42,
            "africa": 0.50,
            "north_america": -1.30,
            "south_america": -3.10,
        },
    }


def test_target_grid_spec_expands_to_the_canonical_global_one_degree_centers():
    """Catch an off-by-half-cell grid shift between products and continent masks."""
    paper = _load("figure03_paper_mascon.json")

    lat, lon = target_grid_coordinates(paper["target_grid"])

    assert lat.shape == (180,)
    assert lon.shape == (360,)
    assert lat[[0, -1]].tolist() == [-89.5, 89.5]
    assert lon[[0, -1]].tolist() == [-179.5, 179.5]
    assert np.diff(lat).tolist() == [1.0] * 179
    assert np.diff(lon).tolist() == [1.0] * 359
