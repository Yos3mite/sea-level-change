from pathlib import Path

import pytest

from pysrc.reference_products.figure3.provenance import (
    build_source_manifest,
    sha256_file,
)


def test_build_source_manifest_records_relative_paths_sizes_and_hashes(
    tmp_path: Path,
):
    """Catch manifests that omit immutable file identity or expose host paths."""
    raw = tmp_path / "raw"
    raw.mkdir()
    first = raw / "a.bin"
    first.write_bytes(b"abc")

    manifest = build_source_manifest(
        {"title": "source", "version": "v2", "url": "https://example.test"},
        [first],
        root=tmp_path,
    )

    assert manifest["source"]["version"] == "v2"
    assert manifest["files"] == [
        {
            "path": "raw/a.bin",
            "bytes": 3,
            "sha256": sha256_file(first),
        }
    ]


def test_build_source_manifest_sorts_files_by_relative_posix_path(tmp_path: Path):
    """Catch platform-dependent or caller-dependent manifest ordering."""
    zulu = tmp_path / "z.bin"
    alpha = tmp_path / "nested" / "a.bin"
    alpha.parent.mkdir()
    zulu.write_bytes(b"z")
    alpha.write_bytes(b"a")

    manifest = build_source_manifest({}, [zulu, alpha], root=tmp_path)

    assert [entry["path"] for entry in manifest["files"]] == [
        "nested/a.bin",
        "z.bin",
    ]


def test_build_source_manifest_rejects_files_outside_root(tmp_path: Path):
    """Catch manifests that silently register files outside the declared source tree."""
    root = tmp_path / "source"
    root.mkdir()
    outside = tmp_path / "outside.bin"
    outside.write_bytes(b"outside")

    with pytest.raises(ValueError, match="outside manifest root"):
        build_source_manifest({}, [outside], root=root)
