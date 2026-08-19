from datetime import date
from hashlib import sha256
from pathlib import Path

import pytest

from gmsl_budget.icgem import DownloadSpec, discover_gfc_downloads, download_gfc, parse_gfc


GFC_TEXT = """begin_of_head
modelname TEST
earth_gravity_constant 3.9860044150e+14
radius 6.3781363000e+06
max_degree 2
norm fully_normalized
time_period_of_data 20200101 - 20200131 (mid: 20200116)
end_of_head
gfc 0 0 1.0 0.0 0.0 0.0
gfc 1 0 0.0 0.0 0.0 0.0
gfc 1 1 1.0e-10 -2.0e-10 0.0 0.0
gfc 2 0 -4.84e-4 0.0 0.0 0.0
gfc 2 1 2.0e-10 3.0e-10 0.0 0.0
gfc 2 2 4.0e-10 5.0e-10 0.0 0.0
"""


def test_parse_gfc_reads_dates_normalization_and_coefficients(tmp_path):
    path = tmp_path / "kfilter_DDK1_GSM-2_2020001-2020031_GRFO_UTCSR_BA01_0603.gfc"
    path.write_text(GFC_TEXT, encoding="ascii")
    epoch = parse_gfc(path, lmax=2)
    assert epoch.start == date(2020, 1, 1)
    assert epoch.end == date(2020, 1, 31)
    assert epoch.midpoint == date(2020, 1, 16)
    assert epoch.c[0, 0] == 0.0
    assert epoch.c[2, 0] == pytest.approx(-4.84e-4)
    assert epoch.s[2, 2] == pytest.approx(5.0e-10)
    assert epoch.normalization == "fully_normalized"


def test_parse_gfc_rejects_non_normalized_coefficients(tmp_path):
    path = tmp_path / "GSM-2_2020001-2020031_GRFO_UTCSR_BA01_0603.gfc"
    path.write_text(GFC_TEXT.replace("fully_normalized", "unnormalized"), encoding="ascii")
    with pytest.raises(ValueError, match="fully_normalized"):
        parse_gfc(path)


def test_discovery_filters_ddk1_and_epoch_midpoint(tmp_path):
    html = tmp_path / "listing.html"
    html.write_text(
        """
        <a href="files/kfilter_DDK1_GSM-2_2018275-2018305_GRFO_UTCSR_BA01_0603.gfc">oct</a>
        <a href="files/kfilter_DDK1_GSM-2_2018286-2018316_GRFO_UTCSR_BA01_0603.gfc">oct-late-duplicate</a>
        <a href="files/kfilter_DDK1_GSM-2_2018275-2018305_GRFO_UTCSR_BB01_0603.gfc">oct-degree96</a>
        <a href="files/kfilter_DDK2_GSM-2_2018306-2018335_GRFO_UTCSR_BA01_0603.gfc">nov-ddk2</a>
        <a href="files/kfilter_DDK1_GSM-2_2018306-2018335_GRFO_UTCSR_BA01_0603.gfc">nov</a>
        """,
        encoding="utf-8",
    )
    specs = discover_gfc_downloads([html.as_uri()], start_month="2018-10", end_month="2018-10")
    assert [spec.filename for spec in specs] == [
        "kfilter_DDK1_GSM-2_2018275-2018305_GRFO_UTCSR_BA01_0603.gfc"
    ]


def test_download_rejects_html_saved_as_gfc(tmp_path):
    source = tmp_path / "error.gfc"
    source.write_text("<html>server error</html>", encoding="utf-8")
    spec = DownloadSpec(source.as_uri(), "error.gfc", "test-page")
    with pytest.raises(ValueError, match="not a valid GFC"):
        download_gfc(spec, tmp_path / "downloads")
    assert not (tmp_path / "downloads" / "error.gfc").exists()


def test_download_validates_expected_sha256_and_writes_atomically(tmp_path):
    source = tmp_path / "source.gfc"
    source.write_text(GFC_TEXT, encoding="ascii")
    expected = sha256(source.read_bytes()).hexdigest()
    spec = DownloadSpec(source.as_uri(), "target.gfc", "test-page", expected_sha256=expected)
    result = download_gfc(spec, tmp_path / "downloads")
    assert result.path.name == "target.gfc"
    assert result.sha256 == expected
    assert result.path.read_bytes() == source.read_bytes()
    assert not Path(str(result.path) + ".part").exists()


def test_download_rejects_hash_mismatch(tmp_path):
    source = tmp_path / "source.gfc"
    source.write_text(GFC_TEXT, encoding="ascii")
    spec = DownloadSpec(source.as_uri(), "target.gfc", "test-page", expected_sha256="0" * 64)
    with pytest.raises(ValueError, match="SHA-256"):
        download_gfc(spec, tmp_path / "downloads")
