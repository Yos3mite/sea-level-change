from __future__ import annotations

from dataclasses import dataclass
from datetime import date, timedelta
from hashlib import sha256
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import quote, urljoin, urlparse, urlsplit, urlunsplit
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd


_EPOCH_PATTERN = re.compile(r"GSM-2_(\d{7})-(\d{7})_")


@dataclass(frozen=True)
class GfcEpoch:
    path: Path
    start: date
    end: date
    midpoint: date
    c: np.ndarray
    s: np.ndarray
    lmax: int
    normalization: str
    radius_m: float
    gravity_constant_m3_s2: float


@dataclass(frozen=True)
class DownloadSpec:
    url: str
    filename: str
    source_page: str
    expected_sha256: str | None = None


@dataclass(frozen=True)
class DownloadResult:
    path: Path
    sha256: str
    size_bytes: int
    source_url: str
    source_page: str


class _LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        if tag.lower() != "a":
            return
        attributes = dict(attrs)
        if "href" in attributes:
            self.links.append(attributes["href"])


def _date_from_year_day(value: str) -> date:
    year = int(value[:4])
    day_of_year = int(value[4:])
    return date(year, 1, 1) + timedelta(days=day_of_year - 1)


def epoch_dates_from_filename(filename: str) -> tuple[date, date, date]:
    match = _EPOCH_PATTERN.search(filename)
    if not match:
        raise ValueError(f"GFC filename does not contain epoch dates: {filename}")
    start = _date_from_year_day(match.group(1))
    end = _date_from_year_day(match.group(2))
    if end < start:
        raise ValueError(f"GFC epoch ends before it starts: {filename}")
    midpoint = start + (end - start) / 2
    return start, end, midpoint


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_gfc(path: str | Path, lmax: int | None = None) -> GfcEpoch:
    source_path = Path(path)
    text = source_path.read_text(encoding="utf-8", errors="strict")
    if "end_of_head" not in text or not any(line.startswith("gfc ") for line in text.splitlines()):
        raise ValueError(f"not a valid GFC file: {source_path}")
    header_text, coefficient_text = text.split("end_of_head", 1)
    header = {}
    for raw_line in header_text.splitlines():
        parts = raw_line.strip().split(maxsplit=1)
        if len(parts) == 2:
            header[parts[0].lower()] = parts[1].strip()
    normalization = header.get("norm", "").lower()
    if normalization != "fully_normalized":
        raise ValueError("GFC coefficients must be fully_normalized")
    available_lmax = int(header.get("max_degree", "0"))
    selected_lmax = available_lmax if lmax is None else int(lmax)
    if selected_lmax < 0 or selected_lmax > available_lmax:
        raise ValueError(f"requested lmax {selected_lmax} exceeds available degree {available_lmax}")
    c = np.zeros((selected_lmax + 1, selected_lmax + 1), dtype=np.float64)
    s = np.zeros_like(c)
    for raw_line in coefficient_text.splitlines():
        parts = raw_line.strip().split()
        if len(parts) < 5 or parts[0] != "gfc":
            continue
        degree, order = int(parts[1]), int(parts[2])
        if degree <= selected_lmax and order <= degree:
            c[degree, order] = _float(parts[3])
            s[degree, order] = _float(parts[4])
    c[0, 0] = 0.0
    c.setflags(write=False)
    s.setflags(write=False)
    start, end, midpoint = epoch_dates_from_filename(source_path.name)
    return GfcEpoch(
        path=source_path.resolve(),
        start=start,
        end=end,
        midpoint=midpoint,
        c=c,
        s=s,
        lmax=selected_lmax,
        normalization=normalization,
        radius_m=_float(header.get("radius", "6378136.3")),
        gravity_constant_m3_s2=_float(header.get("earth_gravity_constant", "3.986004415e14")),
    )


def _read_url(url: str) -> bytes:
    scheme = urlparse(url).scheme.lower()
    if scheme not in {"https", "file"}:
        raise ValueError(f"only HTTPS downloads are allowed: {url}")
    request = Request(url, headers={"User-Agent": "gmsl-budget/0.1"}) if scheme == "https" else url
    with urlopen(request, timeout=60) as response:
        return response.read()


def discover_gfc_downloads(
    series_pages: list[str] | tuple[str, ...],
    start_month: str,
    end_month: str,
) -> list[DownloadSpec]:
    start_period = pd.Period(start_month, freq="M")
    end_period = pd.Period(end_month, freq="M")
    if start_period > end_period:
        raise ValueError("start_month must not be after end_month")
    discovered: dict[pd.Period, tuple[tuple[int, int, str], DownloadSpec]] = {}
    for page in series_pages:
        parser = _LinkParser()
        parser.feed(_read_url(page).decode("utf-8", errors="replace"))
        for href in parser.links:
            filename = Path(urlparse(href).path).name
            if (
                not filename.endswith(".gfc")
                or "kfilter_DDK1_GSM-2_" not in filename
                or "_BA01_" not in filename
            ):
                continue
            try:
                epoch_start, epoch_end, midpoint = epoch_dates_from_filename(filename)
            except ValueError:
                continue
            midpoint_period = pd.Period(midpoint, freq="M")
            if start_period <= midpoint_period <= end_period:
                target = date(midpoint.year, midpoint.month, 15)
                score = (abs((midpoint - target).days), -(epoch_end - epoch_start).days, filename)
                joined = urlsplit(urljoin(page, href))
                encoded_url = urlunsplit(
                    (joined.scheme, joined.netloc, quote(joined.path, safe="/%:@"), joined.query, joined.fragment)
                )
                candidate = DownloadSpec(encoded_url, filename, page)
                if midpoint_period not in discovered or score < discovered[midpoint_period][0]:
                    discovered[midpoint_period] = (score, candidate)
    return [discovered[period][1] for period in sorted(discovered)]


def _validate_gfc_bytes(content: bytes) -> None:
    prefix = content[:256].lstrip().lower()
    if prefix.startswith(b"<html") or b"end_of_head" not in content or b"\ngfc " not in content:
        raise ValueError("downloaded content is not a valid GFC file")


def download_gfc(spec: DownloadSpec, destination: str | Path) -> DownloadResult:
    destination_path = Path(destination)
    destination_path.mkdir(parents=True, exist_ok=True)
    final_path = destination_path / spec.filename
    part_path = destination_path / f"{spec.filename}.part"
    content = _read_url(spec.url)
    _validate_gfc_bytes(content)
    digest = sha256(content).hexdigest()
    if spec.expected_sha256 is not None and digest.lower() != spec.expected_sha256.lower():
        raise ValueError(f"downloaded GFC SHA-256 mismatch for {spec.filename}")
    try:
        part_path.write_bytes(content)
        part_path.replace(final_path)
    finally:
        if part_path.exists():
            part_path.unlink()
    return DownloadResult(final_path.resolve(), digest, len(content), spec.url, spec.source_page)
