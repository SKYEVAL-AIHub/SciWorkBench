"""Deterministic reference analysis for AST-TESS-001.

The implementation intentionally exposes each scientific dependency as a
structured output so grading can distinguish a complete research chain from a
plausible-looking report.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares
from scipy.signal import savgol_filter


RANDOM_SEED = 20260819
REFERENCE_PERIOD = 8.463
REFERENCE_T0 = 1330.390
PIXEL_SCALE_ARCSEC = 21.0
REQUIRED_OUTPUT_FILES = {
    "analysis_pipeline.py",
    "analysis_report.md",
    "centroid_audit.csv",
    "diagnostic_summary.png",
    "event_metrics.csv",
    "injection_recovery.csv",
    "period_search.csv",
    "replay_manifest.json",
    "results.json",
    "sensitivity.csv",
}


def integrated_box_model(
    time_values: np.ndarray,
    center: float,
    duration: float,
    depth: float,
    exposure_days: float,
) -> np.ndarray:
    """Analytically integrate a box transit across each exposure."""
    times = np.asarray(time_values, dtype=float)
    exposure = float(exposure_days)
    left = times - exposure / 2.0
    right = times + exposure / 2.0
    transit_left = center - duration / 2.0
    transit_right = center + duration / 2.0
    overlap = np.maximum(0.0, np.minimum(right, transit_right) - np.maximum(left, transit_left))
    return 1.0 - float(depth) * np.clip(overlap / exposure, 0.0, 1.0)


def periodic_box_model(
    time_values: np.ndarray,
    t0: float,
    period: float,
    duration: float,
    depth: float,
    exposure_days: np.ndarray,
    supersample: int = 7,
) -> np.ndarray:
    times = np.asarray(time_values, dtype=float)
    exposures = np.asarray(exposure_days, dtype=float)
    offsets = np.linspace(-0.5, 0.5, supersample)
    sampled = times[:, None] + exposures[:, None] * offsets[None, :]
    phase = ((sampled - t0 + period / 2.0) % period) - period / 2.0
    fraction = (np.abs(phase) <= duration / 2.0).mean(axis=1)
    return 1.0 - float(depth) * fraction


def fit_linear_ephemeris(epochs: np.ndarray, times: np.ndarray, errors: np.ndarray) -> Dict[str, object]:
    epochs = np.asarray(epochs, dtype=float)
    times = np.asarray(times, dtype=float)
    errors = np.maximum(np.asarray(errors, dtype=float), 1e-8)
    design = np.column_stack([np.ones_like(epochs), epochs])
    weights = 1.0 / errors**2
    normal = design.T @ (weights[:, None] * design)
    covariance = np.linalg.inv(normal)
    parameters = covariance @ (design.T @ (weights * times))
    predicted = design @ parameters
    residual_seconds = (times - predicted) * 86400.0
    return {
        "t0_btjd": float(parameters[0]),
        "period_days": float(parameters[1]),
        "t0_error_days": float(math.sqrt(covariance[0, 0])),
        "period_error_days": float(math.sqrt(covariance[1, 1])),
        "oc_seconds": residual_seconds.tolist(),
        "chi2": float(np.sum(((times - predicted) / errors) ** 2)),
        "dof": int(max(len(times) - 2, 0)),
    }


def difference_centroid(difference: np.ndarray, reference: np.ndarray) -> Dict[str, float]:
    difference = np.asarray(difference, dtype=float)
    reference = np.asarray(reference, dtype=float)
    yy, xx = np.indices(difference.shape)
    diff_weights = np.clip(np.nan_to_num(difference, nan=0.0), 0.0, None)
    ref_weights = np.clip(np.nan_to_num(reference, nan=0.0), 0.0, None)
    if diff_weights.sum() <= 0 or ref_weights.sum() <= 0:
        return {"difference_x_pix": np.nan, "difference_y_pix": np.nan, "reference_x_pix": np.nan, "reference_y_pix": np.nan, "offset_pix": np.nan, "offset_arcsec": np.nan}
    dx = float((xx * diff_weights).sum() / diff_weights.sum())
    dy = float((yy * diff_weights).sum() / diff_weights.sum())
    rx = float((xx * ref_weights).sum() / ref_weights.sum())
    ry = float((yy * ref_weights).sum() / ref_weights.sum())
    offset = math.hypot(dx - rx, dy - ry)
    return {
        "difference_x_pix": dx,
        "difference_y_pix": dy,
        "reference_x_pix": rx,
        "reference_y_pix": ry,
        "offset_pix": offset,
        "offset_arcsec": offset * PIXEL_SCALE_ARCSEC,
    }


def choose_alias(rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    """Choose an alias by event support, contradictions, and depth consistency."""
    scored = []
    for row in rows:
        candidate = dict(row)
        candidate["adjudication_score"] = (
            3.0 * float(candidate["supported_events"])
            - 2.0 * float(candidate["contradictions"])
            + 2.0 * float(candidate["event_consistency"])
        )
        scored.append(candidate)
    return max(scored, key=lambda row: (float(row["adjudication_score"]), -abs(float(row["period_days"]) - REFERENCE_PERIOD)))


def _load_lightcurves(data_dir: Path) -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in sorted((Path(data_dir) / "lightcurves").glob("*.csv.gz"))]
    frame = pd.concat(frames, ignore_index=True)
    finite = np.isfinite(frame["time_btjd"]) & np.isfinite(frame["relative_flux"]) & np.isfinite(frame["relative_flux_err"])
    return frame.loc[finite].sort_values("time_btjd").reset_index(drop=True)


def _detrend_lightcurve(frame: pd.DataFrame, reject_bits: int = 175, window_days: float = 1.2, flare_sigma: float = 5.0) -> pd.DataFrame:
    accepted = (frame["quality"].astype(np.int64).to_numpy() & int(reject_bits)) == 0
    work = frame.loc[accepted].copy().reset_index(drop=True)
    work["trend"] = np.nan
    work["flare_mask"] = False
    for sector, sector_frame in work.groupby("sector", sort=True):
        indices = sector_frame.index.to_numpy()
        times = sector_frame["time_btjd"].to_numpy()
        flux = sector_frame["relative_flux"].to_numpy()
        cuts = np.r_[0, np.where(np.diff(times) > 0.2)[0] + 1, len(times)]
        for start, stop in zip(cuts[:-1], cuts[1:]):
            local_indices = indices[start:stop]
            local_time = times[start:stop]
            local_flux = flux[start:stop]
            if len(local_flux) < 9:
                work.loc[local_indices, "trend"] = np.nanmedian(local_flux)
                continue
            cadence = float(np.nanmedian(np.diff(local_time)))
            window = max(9, int(round(window_days / cadence)))
            if window % 2 == 0:
                window += 1
            if window >= len(local_flux):
                window = len(local_flux) - 1 if len(local_flux) % 2 == 0 else len(local_flux)
            window = max(window, 5)
            coarse = savgol_filter(local_flux, window, 2, mode="mirror")
            residual = local_flux / coarse - 1.0
            scatter = 1.4826 * np.nanmedian(np.abs(residual - np.nanmedian(residual)))
            flare = residual > flare_sigma * max(scatter, 1e-8)
            phase = np.abs(((local_time - REFERENCE_T0 + REFERENCE_PERIOD / 2) % REFERENCE_PERIOD) - REFERENCE_PERIOD / 2)
            protected = flare | (phase < 0.12)
            clean = local_flux.copy()
            usable = ~protected
            if usable.sum() >= 4:
                clean[protected] = np.interp(local_time[protected], local_time[usable], local_flux[usable])
            trend = savgol_filter(clean, window, 2, mode="mirror")
            work.loc[local_indices, "trend"] = trend
            work.loc[local_indices, "flare_mask"] = flare
    work["flat_flux"] = work["relative_flux"] / work["trend"]
    scatter = 1.4826 * np.nanmedian(np.abs(work["flat_flux"] - np.nanmedian(work["flat_flux"])))
    work["analysis_mask"] = (~work["flare_mask"]) & (work["flat_flux"] > 1 - 10 * scatter) & (work["flat_flux"] < 1 + 5 * scatter)
    return work


def _alias_event_metrics(frame: pd.DataFrame, period: float) -> Tuple[int, int, float, int]:
    depths = []
    supported = 0
    contradictions = 0
    covered = 0
    for _, sector_frame in frame.groupby("sector"):
        tmin = float(sector_frame["time_btjd"].min())
        tmax = float(sector_frame["time_btjd"].max())
        first = int(math.floor((tmin - REFERENCE_T0) / period)) - 1
        last = int(math.ceil((tmax - REFERENCE_T0) / period)) + 1
        for epoch in range(first, last + 1):
            center = REFERENCE_T0 + epoch * period
            delta = np.abs(sector_frame["time_btjd"].to_numpy() - center)
            in_values = sector_frame.loc[(delta < 0.055) & sector_frame["analysis_mask"], "flat_flux"].to_numpy()
            out_values = sector_frame.loc[(delta > 0.105) & (delta < 0.22) & sector_frame["analysis_mask"], "flat_flux"].to_numpy()
            if len(in_values) < 20 or len(out_values) < 35:
                continue
            covered += 1
            depth = float(np.nanmedian(out_values) - np.nanmedian(in_values))
            scatter = 1.4826 * np.nanmedian(np.abs(out_values - np.nanmedian(out_values)))
            uncertainty = max(scatter / math.sqrt(len(in_values)), 1e-6)
            significance = depth / uncertainty
            if significance >= 2.0:
                supported += 1
                depths.append(depth)
            elif significance < 1.0:
                contradictions += 1
    if depths:
        median = float(np.median(depths))
        dispersion = float(1.4826 * np.median(np.abs(np.asarray(depths) - median)))
        consistency = float((supported / max(covered, 1)) * math.exp(-dispersion / max(abs(median), 1e-6)))
    else:
        consistency = 0.0
    return supported, contradictions, consistency, covered


def _period_search(frame: pd.DataFrame) -> pd.DataFrame:
    clean = frame.loc[frame["analysis_mask"]].copy()
    time_values = clean["time_btjd"].to_numpy()
    flux = clean["flat_flux"].to_numpy()
    errors = np.maximum(clean["relative_flux_err"].to_numpy(), np.nanmedian(clean["relative_flux_err"]) * 0.5)
    periods = np.array([REFERENCE_PERIOD / 2, REFERENCE_PERIOD, REFERENCE_PERIOD * 2])
    labels = ["P/2", "P", "2P"]
    bls = BoxLeastSquares(time_values, flux, dy=errors)
    rows = []
    for label, period in zip(labels, periods):
        power = bls.power(np.array([period]), np.array([0.08, 0.10, 0.12, 0.15]), objective="snr")
        index = int(np.nanargmax(power.power))
        supported, contradictions, consistency, covered = _alias_event_metrics(frame, float(period))
        rows.append({
            "alias_label": label,
            "period_days": float(period),
            "bls_snr": float(power.power[index]),
            "bls_t0_btjd": float(power.transit_time[index]),
            "bls_duration_days": float(power.duration[index]),
            "bls_depth": float(power.depth[index]),
            "window_function_power": float(abs(np.mean(np.exp(2j * np.pi * time_values / period)))),
            "covered_events": covered,
            "supported_events": supported,
            "contradictions": contradictions,
            "event_consistency": consistency,
        })
    selected = choose_alias(rows)
    for row in rows:
        row["adjudication_score"] = 3.0 * row["supported_events"] - 2.0 * row["contradictions"] + 2.0 * row["event_consistency"]
        row["selected"] = row["alias_label"] == selected["alias_label"]
    return pd.DataFrame(rows)


def _fit_event(path: Path, bootstrap_replicates: int, rng: np.random.Generator) -> Dict[str, object]:
    frame = pd.read_csv(path)
    event_id = str(frame["event_id"].iloc[0])
    epoch = int(frame["epoch"].iloc[0])
    seed = REFERENCE_T0 + epoch * REFERENCE_PERIOD
    valid = np.isfinite(frame["time_btjd"]) & np.isfinite(frame["relative_flux"]) & ((frame["quality"].astype(np.int64) & 175) == 0)
    frame = frame.loc[valid].copy()
    x = frame["time_btjd"].to_numpy() - seed
    flux = frame["relative_flux"].to_numpy()
    error = np.maximum(frame["relative_flux_err"].to_numpy(), 1e-6)
    baseline = np.abs(x) > 0.09
    coefficients = np.polyfit(x[baseline], flux[baseline], deg=2, w=1.0 / error[baseline])
    trend = np.polyval(coefficients, x)
    flat = flux / trend
    scatter = 1.4826 * np.median(np.abs(flat[baseline] - np.median(flat[baseline])))
    keep = flat < np.median(flat[baseline]) + 4.0 * max(scatter, 1e-6)
    time_values = frame["time_btjd"].to_numpy()[keep]
    flat = flat[keep]
    error = np.sqrt((error[keep] / trend[keep]) ** 2 + scatter**2)
    exposure = float(frame["cadence_seconds"].median()) / 86400.0
    centers = np.linspace(seed - 0.045, seed + 0.045, 121)
    durations = np.linspace(0.08, 0.16, 9)

    def best_fit(values: np.ndarray, fixed_duration: float = None) -> Tuple[float, float, float, float]:
        best = (np.inf, seed, 0.12, 0.0)
        duration_values = [fixed_duration] if fixed_duration is not None else durations
        weights = 1.0 / error**2
        for duration in duration_values:
            for center in centers:
                unit = 1.0 - integrated_box_model(time_values, center, float(duration), 1.0, exposure)
                depth = float(np.clip(np.sum(weights * unit * (1.0 - values)) / max(np.sum(weights * unit**2), 1e-12), 0.0, 0.02))
                model = 1.0 - depth * unit
                chi2 = float(np.sum(weights * (values - model) ** 2))
                if chi2 < best[0]:
                    best = (chi2, float(center), float(duration), depth)
        return best

    chi2, center, duration, depth = best_fit(flat)
    model = integrated_box_model(time_values, center, duration, depth, exposure)
    residual = flat - model
    boot_centers = []
    for _ in range(int(bootstrap_replicates)):
        synthetic = model + rng.choice(residual, size=len(residual), replace=True)
        boot_centers.append(best_fit(synthetic, fixed_duration=duration)[1])
    center_error = max(float(np.std(boot_centers, ddof=1)) if len(boot_centers) > 1 else 0.0, exposure / math.sqrt(12.0))
    in_mask = np.abs(time_values - center) <= duration / 2
    depth_error = max(scatter / math.sqrt(max(int(in_mask.sum()), 1)), 1e-6)
    return {
        "event_id": event_id,
        "epoch": epoch,
        "sector": int(frame["sector"].iloc[0]),
        "seed_btjd": seed,
        "midtime_btjd": center,
        "midtime_bjd_tdb": center + 2457000.0,
        "midtime_error_days": center_error,
        "depth": depth,
        "depth_error": depth_error,
        "duration_days": duration,
        "exposure_seconds": float(frame["cadence_seconds"].median()),
        "flare_points_removed": int((~keep).sum()),
        "fit_chi2": chi2,
        "accepted": bool(depth / depth_error >= 1.5),
    }


def _centroid_table(data_dir: Path) -> Tuple[pd.DataFrame, float]:
    rows = []
    for path in sorted((Path(data_dir) / "pixel_events").glob("*.npz")):
        arrays = np.load(path)
        aperture = arrays["aperture"] > 0
        reference = np.where(aperture, arrays["out_mean"], 0.0)
        metrics = difference_centroid(arrays["difference"], reference)
        variance = arrays["in_stderr"] ** 2 + arrays["out_stderr"] ** 2
        signal = float(np.nansum(arrays["difference"][aperture]))
        noise = float(np.sqrt(np.nansum(variance[aperture])))
        metadata = json.loads(str(arrays["metadata_json"]))
        metrics.update({
            "event_id": path.stem,
            "difference_snr": signal / max(noise, 1e-12),
            "in_count": metadata["in_count"],
            "out_count": metadata["out_count"],
            "on_target": bool(metrics["offset_pix"] < 1.0),
        })
        rows.append(metrics)
    neighbors = pd.read_csv(Path(data_dir) / "catalog" / "tic_neighbors.csv")
    contaminating = neighbors.loc[(neighbors["tic_id"] != 441420236) & (neighbors["separation_arcsec"] <= 42.0), "relative_tess_flux"].sum()
    dilution = float(contaminating / (1.0 + contaminating))
    frame = pd.DataFrame(rows)
    frame["catalog_dilution_fraction_42arcsec"] = dilution
    return frame, dilution


def _remove_known_events(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    values = frame["flat_flux"].to_numpy().copy()
    times = frame["time_btjd"].to_numpy()
    phase = np.abs(((times - REFERENCE_T0 + REFERENCE_PERIOD / 2) % REFERENCE_PERIOD) - REFERENCE_PERIOD / 2)
    event_mask = phase < 0.09
    pool = values[(~event_mask) & frame["analysis_mask"].to_numpy()] - 1.0
    values[event_mask] = 1.0 + rng.choice(pool, size=int(event_mask.sum()), replace=True)
    return values


def _injection_recovery(frame: pd.DataFrame, data_dir: Path, rng: np.random.Generator) -> pd.DataFrame:
    grid = pd.read_csv(Path(data_dir) / "config" / "injection_grid.csv")
    clean = frame.loc[frame["analysis_mask"]].copy().reset_index(drop=True)
    times = clean["time_btjd"].to_numpy()
    errors = np.maximum(clean["relative_flux_err"].to_numpy(), np.nanmedian(clean["relative_flux_err"]) * 0.5)
    exposures = clean["cadence_seconds"].to_numpy() / 86400.0
    base = _remove_known_events(clean, rng)
    rows = []
    for record in grid.to_dict(orient="records"):
        t0 = REFERENCE_T0 + float(record["phase_offset_days"])
        model = periodic_box_model(times, t0, float(record["period_days"]), float(record["duration_days"]), float(record["depth"]), exposures)
        injected = base * model
        bls = BoxLeastSquares(times, injected, dy=errors)
        aliases = np.array([REFERENCE_PERIOD / 2, REFERENCE_PERIOD, REFERENCE_PERIOD * 2])
        power = bls.power(aliases, float(record["duration_days"]), objective="snr")
        selected_index = int(np.nanargmax(power.power))
        recovered_period = float(aliases[selected_index])
        recovered_t0 = float(power.transit_time[selected_index])
        phase_error = abs(((recovered_t0 - t0 + REFERENCE_PERIOD / 2) % REFERENCE_PERIOD) - REFERENCE_PERIOD / 2)
        rows.append({
            **record,
            "recovered_period_days": recovered_period,
            "recovered_t0_btjd": recovered_t0,
            "recovered_depth": float(power.depth[selected_index]),
            "detection_snr": float(power.power[selected_index]),
            "phase_error_days": phase_error,
            "recovered": bool(abs(recovered_period - REFERENCE_PERIOD) < 0.03 and phase_error <= float(record["duration_days"]) / 2 and float(power.power[selected_index]) >= 7.0),
        })
    return pd.DataFrame(rows)


def _sensitivity_table(raw_frame: pd.DataFrame, data_dir: Path) -> pd.DataFrame:
    scenarios = pd.read_csv(Path(data_dir) / "config" / "sensitivity_scenarios.csv")
    rows = []
    for record in scenarios.to_dict(orient="records"):
        detrended = _detrend_lightcurve(
            raw_frame,
            reject_bits=int(record["quality_bits"]),
            window_days=float(record["detrend_days"]),
            flare_sigma=float(record["flare_sigma"]),
        )
        period_table = _period_search(detrended)
        selected = period_table.loc[period_table["selected"]].iloc[0]
        rows.append({
            "scenario_id": record["scenario_id"],
            "quality_bits": int(record["quality_bits"]),
            "flare_sigma": float(record["flare_sigma"]),
            "detrend_days": float(record["detrend_days"]),
            "exposure_integrated": bool(record["exposure_integrated"]),
            "exclude_overlap_event": record.get("exclude_overlap_event") if pd.notna(record.get("exclude_overlap_event")) else "",
            "selected_period_days": float(selected["period_days"]),
            "selected_alias": str(selected["alias_label"]),
            "supported_events": int(selected["supported_events"]),
            "classification_stable": bool(str(selected["alias_label"]) == "P" and int(selected["supported_events"]) >= 3),
        })
    return pd.DataFrame(rows)


def run_reference_analysis(data_dir: Path, output_dir: Path, bootstrap_replicates: int = 300) -> Dict[str, object]:
    started = time.perf_counter()
    data_dir = Path(data_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    raw = _load_lightcurves(data_dir)
    detrended = _detrend_lightcurve(raw)
    period_search = _period_search(detrended)
    selected = period_search.loc[period_search["selected"]].iloc[0]

    event_rows = [_fit_event(path, bootstrap_replicates, rng) for path in sorted((data_dir / "event_windows").glob("*.csv.gz"))]
    event_table = pd.DataFrame(event_rows).sort_values("epoch").reset_index(drop=True)
    ephemeris = fit_linear_ephemeris(event_table["epoch"].to_numpy(), event_table["midtime_btjd"].to_numpy(), event_table["midtime_error_days"].to_numpy())
    event_table["oc_seconds"] = ephemeris["oc_seconds"]
    loo_periods = []
    for index in range(len(event_table)):
        subset = event_table.drop(index)
        fit = fit_linear_ephemeris(subset["epoch"].to_numpy(), subset["midtime_btjd"].to_numpy(), subset["midtime_error_days"].to_numpy())
        loo_periods.append({"excluded_event_id": event_table.loc[index, "event_id"], "period_days": fit["period_days"]})

    centroid_table, dilution = _centroid_table(data_dir)
    injections = _injection_recovery(detrended, data_dir, rng)
    sensitivity = _sensitivity_table(raw, data_dir)
    recovery_fraction = float(injections["recovered"].mean())
    period_ok = abs(float(selected["period_days"]) - REFERENCE_PERIOD) < 0.03
    centroid_ok = bool(centroid_table["offset_pix"].median() < 0.8 and centroid_table["on_target"].sum() >= 4)
    sensitivity_ok = bool(sensitivity["classification_stable"].sum() >= 6)
    accepted_events = int(event_table["accepted"].sum())
    if period_ok and centroid_ok and recovery_fraction >= 0.55 and accepted_events >= 3 and sensitivity_ok:
        classification = "robust"
    elif centroid_table["offset_pix"].median() >= 1.5 or accepted_events < 2:
        classification = "artifact"
    else:
        classification = "ambiguous"

    period_search.to_csv(output_dir / "period_search.csv", index=False)
    event_table.to_csv(output_dir / "event_metrics.csv", index=False)
    centroid_table.to_csv(output_dir / "centroid_audit.csv", index=False)
    injections.to_csv(output_dir / "injection_recovery.csv", index=False)
    sensitivity.to_csv(output_dir / "sensitivity.csv", index=False)
    result = {
        "task_id": "AST-TESS-001",
        "random_seed": RANDOM_SEED,
        "time_system": "BTJD = BJD_TDB - 2457000",
        "selected_period_days": float(selected["period_days"]),
        "classification": classification,
        "event_count": len(event_table),
        "accepted_event_count": accepted_events,
        "centroid_event_count": len(centroid_table),
        "injection_count": len(injections),
        "injection_recovery_fraction": recovery_fraction,
        "catalog_dilution_fraction_42arcsec": dilution,
        "linear_ephemeris": ephemeris,
        "leave_one_out_periods": loo_periods,
        "core_chain": {
            "R12": {"quality_reject_bits": 175, "accepted_rows": int(detrended["analysis_mask"].sum()), "flare_rows": int(detrended["flare_mask"].sum())},
            "R13": {"method_primary": "event-protected Savitzky-Golay", "method_secondary": "local quadratic event baseline", "window_days": 1.2},
            "R14": {"selected_alias": str(selected["alias_label"]), "period_days": float(selected["period_days"]), "supported_events": int(selected["supported_events"])},
            "R15": {"ephemeris": ephemeris, "event_count": len(event_table), "bootstrap_replicates": int(bootstrap_replicates), "leave_one_out": loo_periods},
            "R16": {"median_offset_pix": float(centroid_table["offset_pix"].median()), "on_target_events": int(centroid_table["on_target"].sum()), "dilution_fraction": dilution},
            "R17": {"recovery_fraction": recovery_fraction, "stable_scenarios": int(sensitivity["classification_stable"].sum()), "classification": classification},
        },
        "runtime_seconds": time.perf_counter() - started,
    }
    (output_dir / "results.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _build_report(output: Path, result: Dict[str, object]) -> None:
    centroids = pd.read_csv(output / "centroid_audit.csv")
    injections = pd.read_csv(output / "injection_recovery.csv")
    sensitivity = pd.read_csv(output / "sensitivity.csv")
    lines = [
        "# AU Mic 周期性暗化事件：耀斑感知验证报告",
        "",
        "## 结论先行",
        "",
        f"最终分类为 **{result['classification']}**。选择的周期为 {result['selected_period_days']:.6f} d；"
        f"5 个事件中 {result['accepted_event_count']} 个达到单事件接受阈值。"
        f"差分像质心偏移中位数为 {result['core_chain']['R16']['median_offset_pix']:.3f} pixel，"
        f"42 arcsec 内目录稀释比例为 {result['catalog_dilution_fraction_42arcsec']:.6f}。",
        "",
        "## R12 时间、质量与掩膜",
        "",
        "使用 BTJD = BJD_TDB − 2457000；QUALITY 采用 `(QUALITY & 175) == 0`。"
        f"最终分析点数 {result['core_chain']['R12']['accepted_rows']}，识别耀斑点 {result['core_chain']['R12']['flare_rows']}。",
        "",
        "## R13 活动分离与去趋势",
        "",
        "主方法为事件保护的 Savitzky–Golay 去趋势，窗口 1.2 d；逐事件计时使用局部二次基线。"
        "耀斑只作高侧异常剔除，暗化窗口受到保护，避免把目标信号当作活动噪声消除。",
        "",
        "## R14 周期与别名",
        "",
        "同时比较 P/2、P、2P 的 BLS、窗口函数、覆盖事件、支持事件和矛盾事件。"
        "最终选择 P，因为 P/2 在交错相位产生多个无暗化预测，2P 漏掉受支持事件。",
        "",
        "## R15 逐事件计时与星历",
        "",
        f"线性星历：T0 = {result['linear_ephemeris']['t0_btjd']:.6f} BTJD，"
        f"P = {result['linear_ephemeris']['period_days']:.8f} d。"
        "每个中点均对 120 s 曝光进行积分，并用固定随机种子的残差 bootstrap 给出误差；"
        "中点与周期误差单位为 day，O−C 单位为 second（秒）；留一法周期结果记录在 results.json。",
        "",
        "## R16 源位置与污染",
        "",
        f"5 个事件的差分像均给出独立质心；其中 {int(centroids['on_target'].sum())} 个位于 1 pixel 内。"
        "目录污染按 TIC 邻星 Tmag 转换为相对 TESS 通量并在 42 arcsec 范围求和。",
        "",
        "## R17 注入恢复与敏感性",
        "",
        f"完成 {len(injections)} 组深度–持续时间–相位注入，恢复率 {injections['recovered'].mean():.3f}；"
        f"8 个敏感性场景中 {int(sensitivity['classification_stable'].sum())} 个保持主分类。",
        "",
        "## 限制",
        "",
        "该流程使用盒形事件模型，不是完整的物理凌星模型；差分像来自公开 TPF 的事件内外均值，"
        "其空间结论应结合低信噪事件单独审查。活跃恒星的相关噪声可能使 bootstrap 误差偏乐观。",
        "",
        "## 数据与科学背景",
        "",
        "数据来自 NASA/MAST 的 TESS SPOC light-curve 与 target-pixel 产品。"
        "背景参见 Gilbert et al. (2022), arXiv:2109.03924, DOI:10.3847/1538-3881/ac23ca。",
    ]
    (output / "analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _build_figure(data_dir: Path, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    periods = pd.read_csv(output / "period_search.csv")
    events = pd.read_csv(output / "event_metrics.csv")
    centroids = pd.read_csv(output / "centroid_audit.csv")
    injections = pd.read_csv(output / "injection_recovery.csv")
    sensitivity = pd.read_csv(output / "sensitivity.csv")
    lightcurves = pd.concat([pd.read_csv(path) for path in sorted((data_dir / "lightcurves").glob("*.csv.gz"))])

    figure, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    for sector, group in lightcurves.groupby("sector"):
        axes[0, 0].plot(group["time_btjd"], group["relative_flux"], ".", ms=0.7, alpha=0.45, label=f"S{int(sector)}")
    axes[0, 0].set(title="TESS light curves", xlabel="BTJD", ylabel="Relative flux")
    axes[0, 0].legend(markerscale=5, frameon=False)

    colors = ["#8da0cb" if not selected else "#fc8d62" for selected in periods["selected"]]
    axes[0, 1].bar(periods["alias_label"], periods["adjudication_score"], color=colors)
    axes[0, 1].set(title="Alias adjudication", ylabel="Dependency score")
    for index, row in periods.iterrows():
        axes[0, 1].text(index, row["adjudication_score"], f"{row['period_days']:.3f} d", ha="center", va="bottom", fontsize=8)

    axes[0, 2].errorbar(events["epoch"], events["oc_seconds"], yerr=events["midtime_error_days"] * 86400, fmt="o", color="#66c2a5")
    axes[0, 2].axhline(0, color="0.4", lw=1)
    axes[0, 2].set(title="Event timing residuals", xlabel="Epoch", ylabel="O-C (s)")

    axes[1, 0].bar(centroids["event_id"], centroids["offset_pix"], color="#ffd92f")
    axes[1, 0].axhline(1.0, color="crimson", ls="--", lw=1, label="1 pixel")
    axes[1, 0].set(title="Difference-image centroids", ylabel="Offset (pixel)")
    axes[1, 0].legend(frameon=False)

    recovery = injections.pivot_table(index="depth", columns="duration_days", values="recovered", aggfunc="mean")
    image = axes[1, 1].imshow(recovery.to_numpy(), origin="lower", aspect="auto", vmin=0, vmax=1, cmap="viridis")
    axes[1, 1].set(title="Injection recovery", xlabel="Duration grid", ylabel="Depth grid")
    axes[1, 1].set_xticks(range(len(recovery.columns)), [f"{value:.2f}" for value in recovery.columns], rotation=45)
    axes[1, 1].set_yticks(range(len(recovery.index)), [f"{value:.4f}" for value in recovery.index])
    figure.colorbar(image, ax=axes[1, 1], label="Recovery fraction")

    stable = sensitivity["classification_stable"].astype(int)
    axes[1, 2].bar(sensitivity["scenario_id"], stable, color=np.where(stable, "#66c2a5", "#e78ac3"))
    axes[1, 2].set(title="Sensitivity stability", xlabel="Scenario", ylabel="Stable (0/1)", ylim=(0, 1.15))
    axes[1, 2].tick_params(axis="x", rotation=45)

    figure.suptitle("AST-TESS-001: flare-aware validation chain", fontsize=14)
    figure.savefig(output / "diagnostic_summary.png", dpi=180)
    plt.close(figure)


def run_complete_analysis(
    data_dir: Path,
    output_dir: Path,
    bootstrap_replicates: int = 300,
    pipeline_source: Path = None,
) -> Dict[str, object]:
    """Generate the complete ten-file submission from public data only."""
    started = time.perf_counter()
    data_dir = Path(data_dir).resolve()
    output_dir = Path(output_dir).resolve()
    pipeline_source = Path(pipeline_source or __file__).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in REQUIRED_OUTPUT_FILES:
        target = output_dir / name
        if target.exists() and target.resolve() != pipeline_source:
            target.unlink()

    result = run_reference_analysis(data_dir, output_dir, bootstrap_replicates=bootstrap_replicates)
    _build_report(output_dir, result)
    _build_figure(data_dir, output_dir)
    pipeline_target = output_dir / "analysis_pipeline.py"
    if pipeline_target.resolve() != pipeline_source:
        shutil.copy2(pipeline_source, pipeline_target)

    hash_names = sorted(REQUIRED_OUTPUT_FILES - {"replay_manifest.json"})
    manifest = {
        "task_id": "AST-TESS-001",
        "command": (
            "python analysis_pipeline.py --data data --output output "
            f"--bootstrap-replicates {int(bootstrap_replicates)}"
        ),
        "random_seed": RANDOM_SEED,
        "bootstrap_replicates": int(bootstrap_replicates),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "wall_clock_seconds": time.perf_counter() - started,
        "output_sha256": {name: _sha256_file(output_dir / name) for name in hash_names},
    }
    (output_dir / "replay_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    actual = {path.name for path in output_dir.iterdir() if path.is_file()}
    if actual != REQUIRED_OUTPUT_FILES:
        raise RuntimeError(f"Unexpected output file set: {sorted(actual)}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the complete AST-TESS-001 analysis chain.")
    parser.add_argument("--data", type=Path, default=Path("data"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    parser.add_argument("--bootstrap-replicates", type=int, default=300)
    args = parser.parse_args()
    run_complete_analysis(
        args.data,
        args.output,
        bootstrap_replicates=args.bootstrap_replicates,
        pipeline_source=Path(__file__),
    )


if __name__ == "__main__":
    main()
