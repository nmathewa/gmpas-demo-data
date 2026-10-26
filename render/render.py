"""Render the fetched fields to value tiles for the gmpas demo.

    python render/render.py RAW_DIR MESH_FILE [SITE_DIR]

RAW_DIR holds what fetch.py wrote (RAW_DIR/<field>/<time>.npy). For each zoom
level z the global lon/lat grid is 512*2^z by 256*2^z pixels, and every pixel
takes the value of the mesh cell nearest its centre -- gmpas's own ViewIndex,
the same lookup `gmpas view` draws with. Tiles are 256-pixel PNGs in Leaflet's
EPSG:4326 layout, tile (0, 0) at the north-west corner.

A tile holds values, not colours, so the page can colour them with any
colormap and range, and answer a click with the number:

    byte 0       no value (NaN)
    byte k>0     lo + (k - 1) * (hi - lo) / 254

with lo, hi the field's own minimum and maximum over every rendered step,
written to manifest.json. That is 1/254 of the field's range: about 1 W m^-2
for OLR, 0.4 K for 2 m temperature, 41 Pa for sea-level pressure.

Writes SITE_DIR/tiles/<field>/<time>/<z>/<x>/<y>.png and SITE_DIR/manifest.json.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image

from gmpas.mesh import MpasMesh
from gmpas.viewer import ViewIndex

try:                                     # gmpas > 0.9.0
    from gmpas.data import pretty_units
except ImportError:
    def pretty_units(units: str) -> str:
        return units

HERE = Path(__file__).resolve().parent
TILE = 256
LEVELS = 254                             # byte 0 is "no value"
DATA_URL = "https://data.gdex.ucar.edu/d010129/dyamond2/3.75km"
HOME = (-180.0, 180.0, -90.0, 90.0)
PROJECTIONS = ("auto", "Orthographic", "NorthPolarStereo", "SouthPolarStereo",
               "LambertConformal", "Mercator")

_index: dict[int, ViewIndex] = {}        # zoom -> view index, shared by fork


def attrs(raw: Path, fields, stamp: str) -> dict:
    """long_name and units per field, read once from the remote file."""
    cache = raw / "attrs.json"
    if cache.exists():
        return json.loads(cache.read_text())
    import netCDF4

    with netCDF4.Dataset(f"{DATA_URL}/diag.{stamp}.nc#mode=bytes") as nc:
        out = {f: {"long_name": getattr(nc[f], "long_name", f),
                   "units": getattr(nc[f], "units", "")} for f in fields}
    cache.write_text(json.dumps(out, indent=1))
    return out


def ranges(raw: Path, fields, times) -> dict:
    """Each field's min and max over every step (the byte scale), and its
    2nd/98th percentile (the page's default colour range), cached."""
    cache = raw / "ranges.json"
    out = json.loads(cache.read_text()) if cache.exists() else {}
    for f in fields:
        if f in out:
            continue
        lo = hi = None
        sample = []
        for s in times:
            v = np.load(raw / f / f"{s}.npy")
            v = v[np.isfinite(v)]
            lo = float(v.min()) if lo is None else min(lo, float(v.min()))
            hi = float(v.max()) if hi is None else max(hi, float(v.max()))
            sample.append(v[:: max(1, v.size // 200_000)])
        p2, p98 = np.percentile(np.concatenate(sample), [2, 98])
        out[f] = {"lo": lo, "hi": hi, "p2": float(p2), "p98": float(p98)}
        cache.write_text(json.dumps(out, indent=1))
    return out


def frame(job) -> str:
    raw, site, f, stamp, lo, hi = job
    out = site / "tiles" / f["name"] / stamp
    if (out / "0" / "0" / "0.png").exists():
        return f"{f['name']} {stamp} cached"
    t = time.time()
    values = np.load(raw / f["name"] / f"{stamp}.npy")
    scale = LEVELS / (hi - lo) if hi > lo else 0.0
    n = 0
    for z in range(f["maxzoom"], -1, -1):
        img = _index[z].frame(values)                      # nearest cell, NaN off
        code = np.where(np.isfinite(img),
                        1 + np.clip(np.rint((img - lo) * scale), 0, LEVELS), 0).astype(np.uint8)
        h, w = code.shape
        for x in range(w // TILE):
            for y in range(h // TILE):
                d = out / str(z) / str(x)
                d.mkdir(parents=True, exist_ok=True)
                Image.fromarray(code[y * TILE:(y + 1) * TILE, x * TILE:(x + 1) * TILE], "L") \
                    .save(d / f"{y}.png", "PNG", compress_level=9)
                n += 1
    return f"{f['name']} {stamp} {n} tiles {time.time() - t:.1f}s"


def homes() -> dict:
    """The globe framed in each projection the demo offers: what the page's
    projection select asks api/projview for."""
    from gmpas.projection import home_view

    out = {}
    for name in PROJECTIONS:
        proj, params, extent = home_view(name, None, HOME, 100.0)
        out[name] = {"proj": proj, "extent": [float(v) for v in extent],
                     "plon": params.get("central_longitude"),
                     "plat": params.get("central_latitude"),
                     "psp": list(params["standard_parallels"])
                     if "standard_parallels" in params else None}
    return out


def main(raw: Path, mesh_file: Path, site: Path) -> None:
    cfg = json.loads((HERE / "fields.json").read_text())
    fields, times = cfg["fields"], cfg["times"]
    names = [f["name"] for f in fields]
    meta = attrs(raw, names, times[0])
    rng = ranges(raw, names, times)

    t = time.time()
    mesh = MpasMesh.load(mesh_file)
    for z in range(max(f["maxzoom"] for f in fields) + 1):
        _index[z] = ViewIndex(mesh, HOME, 512 * 2 ** z, 256 * 2 ** z)
    print(f"indexes for {mesh.n_cells:,} cells in {time.time() - t:.0f}s", flush=True)

    jobs = [(raw, site, f, s, rng[f["name"]]["lo"], rng[f["name"]]["hi"])
            for f in fields for s in times if (raw / f["name"] / f"{s}.npy").exists()]
    with mp.get_context("fork").Pool(8) as pool:
        for line in pool.imap_unordered(frame, jobs):
            print(line, flush=True)

    manifest = {
        "run": cfg["run"],
        "cells": mesh.n_cells,
        "home": list(HOME),
        "coverage": 100.0,
        "times": times,
        "labels": [s.replace("_", " ").replace(".00.00", ":00") for s in times],
        "proj": homes(),
        "tile": {"size": TILE, "format": "png", "levels": LEVELS},
        "fields": [{
            "name": f["name"],
            "label": f"{meta[f['name']]['long_name']} [{pretty_units(meta[f['name']]['units'])}]",
            "lo": rng[f["name"]]["lo"], "hi": rng[f["name"]]["hi"],
            "show": [rng[f["name"]]["p2"], rng[f["name"]]["p98"]],
            "cmap": f.get("cmap", "viridis"),
            "maxzoom": f["maxzoom"],
            "times": [s for s in times
                      if (site / "tiles" / f["name"] / s / "0" / "0" / "0.png").exists()],
        } for f in fields],
    }
    (site / "manifest.json").write_text(json.dumps(manifest, indent=1))
    total = sum(p.stat().st_size for p in (site / "tiles").rglob("*.png"))
    print(f"manifest written; tiles total {total / 1e6:.0f} MB", flush=True)


if __name__ == "__main__":
    site = Path(sys.argv[3]) if len(sys.argv) > 3 else HERE.parent / "site"
    main(Path(sys.argv[1]), Path(sys.argv[2]), site)
