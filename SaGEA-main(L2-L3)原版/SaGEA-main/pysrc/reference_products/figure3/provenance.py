"""Immutable provenance helpers for external Figure 3 inputs."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Iterable, Mapping


_HASH_BLOCK_BYTES = 1024 * 1024


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of *path* using bounded memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(_HASH_BLOCK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def build_source_manifest(
    source: Mapping[str, Any],
    files: Iterable[Path],
    root: Path,
) -> dict[str, Any]:
    """Inventory source files beneath *root* with stable relative paths."""
    root_resolved = Path(root).resolve(strict=True)
    inventory: list[dict[str, Any]] = []

    for file_path in files:
        resolved = Path(file_path).resolve(strict=True)
        if not resolved.is_file():
            raise ValueError(f"manifest entry is not a file: {file_path}")
        try:
            relative = resolved.relative_to(root_resolved)
        except ValueError as exc:
            raise ValueError(f"file is outside manifest root: {file_path}") from exc
        inventory.append(
            {
                "path": relative.as_posix(),
                "bytes": resolved.stat().st_size,
                "sha256": sha256_file(resolved),
            }
        )

    inventory.sort(key=lambda entry: entry["path"])
    return {"source": dict(source), "files": inventory}
