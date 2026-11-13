"""Independent numerical Oracle derived only from public data and submitted tables."""

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_hashes_valid(data_dir):
    data_dir = Path(data_dir)
    try:
        manifest = pd.read_csv(data_dir / "materials_manifest.csv")
    except Exception:
        return False
    if len(manifest) != 55 or not manifest["material_id"].is_unique:
        return False
    for row in manifest.itertuples(index=False):
        path = data_dir / row.relative_path
        if not path.is_file() or sha256_file(path) != row.sha256:
            return False
    return True


def weighted_ephemeris(events):
    epochs = events["epoch"].to_numpy(dtype=float)
    times = events["midtime_btjd"].to_numpy(dtype=float)
    errors = np.maximum(events["midtime_error_days"].to_numpy(dtype=float), 1e-8)
    design = np.column_stack([np.ones(len(events)), epochs])
    weights = 1.0 / errors**2
    covariance = np.linalg.inv(design.T @ (weights[:, None] * design))
    parameters = covariance @ (design.T @ (weights * times))
    predicted = design @ parameters
    return {
        "t0_btjd": float(parameters[0]),
        "period_days": float(parameters[1]),
        "t0_error_days": float(math.sqrt(covariance[0, 0])),
        "period_error_days": float(math.sqrt(covariance[1, 1])),
        "oc_seconds": ((times - predicted) * 86400.0).tolist(),
    }


def leave_one_out(events):
    rows = []
    for index, event in events.reset_index(drop=True).iterrows():
        fit = weighted_ephemeris(events.reset_index(drop=True).drop(index))
        rows.append({"excluded_event_id": event["event_id"], "period_days": fit["period_days"]})
    return rows


def _centroid(difference, reference):
    yy, xx = np.indices(difference.shape)
    dw = np.clip(np.nan_to_num(difference, nan=0.0), 0.0, None)
    rw = np.clip(np.nan_to_num(reference, nan=0.0), 0.0, None)
    dx = float((xx * dw).sum() / dw.sum())
    dy = float((yy * dw).sum() / dw.sum())
    rx = float((xx * rw).sum() / rw.sum())
    ry = float((yy * rw).sum() / rw.sum())
    return dx, dy, rx, ry, math.hypot(dx - rx, dy - ry)


def centroid_expectations(data_dir):
    results = {}
    for path in sorted((Path(data_dir) / "pixel_events").glob("*.npz")):
        arrays = np.load(path)
        aperture = arrays["aperture"] > 0
        reference = np.where(aperture, arrays["out_mean"], 0.0)
        dx, dy, rx, ry, offset = _centroid(arrays["difference"], reference)
        results[path.stem] = {
            "difference_x_pix": dx,
            "difference_y_pix": dy,
            "reference_x_pix": rx,
            "reference_y_pix": ry,
            "offset_pix": offset,
            "offset_arcsec": offset * 21.0,
        }
    return results


def dilution_expectation(data_dir):
    neighbors = pd.read_csv(Path(data_dir) / "catalog" / "tic_neighbors.csv")
    flux = neighbors.loc[(neighbors["tic_id"] != 441420236) & (neighbors["separation_arcsec"] <= 42.0), "relative_tess_flux"].sum()
    return float(flux / (1.0 + flux))


def bool_series(series):
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().map({"true": True, "false": False, "1": True, "0": False})
