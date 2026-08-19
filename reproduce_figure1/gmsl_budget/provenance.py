from __future__ import annotations

from hashlib import sha256
import importlib.metadata
import json
from pathlib import Path
import platform
from typing import Any, Iterable

import numpy as np

from .models import SpatialMask


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_mask(mask: SpatialMask) -> str:
    digest = sha256()
    for value in (mask.latitude, mask.longitude, mask.ocean_fraction):
        array = np.ascontiguousarray(value, dtype="<f8")
        digest.update(str(array.shape).encode("ascii"))
        digest.update(array.tobytes())
    support = np.ascontiguousarray(mask.support, dtype=np.uint8)
    digest.update(str(support.shape).encode("ascii"))
    digest.update(support.tobytes())
    return digest.hexdigest()


def software_versions(packages: Iterable[str] = ()) -> dict[str, str]:
    names = tuple(packages) or ("numpy", "pandas", "scipy", "statsmodels", "xarray", "matplotlib", "h5py")
    result = {"python": platform.python_version(), "platform": platform.platform()}
    for name in names:
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = "not installed"
    return result


def write_json(path: str | Path, payload: Any) -> None:
    Path(path).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
