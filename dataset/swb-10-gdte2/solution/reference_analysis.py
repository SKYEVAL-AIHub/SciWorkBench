"""Deterministic reference pipeline for the GdTe2 boundary-order benchmark.

The script reads only ../data by default. It does not require the paper,
supporting information, original Nanonis files, network access, or a GPU.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import ndimage
from scipy.signal import savgol_filter


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = SCRIPT_DIR.parent / "data"


def robust_limits(image: np.ndarray, percentiles=(1.0, 99.0)) -> tuple[float, float]:
    values = np.nanpercentile(image, percentiles)
    return float(values[0]), float(values[1])


def flatten_linewise(image: np.ndarray) -> np.ndarray:
    x = np.linspace(-1.0, 1.0, image.shape[1])
    design = np.column_stack([x, np.ones_like(x)])
    coefficients = np.linalg.lstsq(design, image.T, rcond=None)[0]
    background = (design @ coefficients).T
    flattened = image - background
    flattened -= np.nanmedian(flattened)
    return flattened


def register_retrace(forward: np.ndarray, retrace: np.ndarray) -> tuple[np.ndarray, list[int], float]:
    a = forward.astype(float) - np.nanmedian(forward)
    b = retrace.astype(float) - np.nanmedian(retrace)
    cross = np.fft.fft2(a) * np.conj(np.fft.fft2(b))
    cross /= np.maximum(np.abs(cross), 1e-30)
    correlation = np.fft.ifft2(cross).real
    shift = np.array(np.unravel_index(np.argmax(correlation), correlation.shape), dtype=int)
    shape = np.array(correlation.shape)
    shift[shift > shape // 2] -= shape[shift > shape // 2]
    registered = np.roll(retrace, shift=tuple(shift), axis=(0, 1)).astype(float)
    corrcoef = float(np.corrcoef(a.ravel(), (registered - np.nanmedian(registered)).ravel())[0, 1])
    return registered, [int(v) for v in shift], corrcoef


def load_frame(path: Path) -> dict:
    with np.load(path) as npz:
        result = {key: npz[key].copy() for key in npz.files}
    registered, shift, correlation = register_retrace(result["forward"], result["retrace"])
    forward = result["forward"].astype(float)
    forward -= np.nanmedian(forward)
    registered -= np.nanmedian(registered)
    result["registered_retrace"] = registered
    result["mean"] = 0.5 * (forward + registered)
    result["registration_shift_pixels_yx"] = shift
    result["forward_retrace_correlation"] = correlation
    return result


def windowed_fft(image: np.ndarray, dx_nm: float, dy_nm: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = flatten_linewise(image)
    wy = np.hanning(data.shape[0])
    wx = np.hanning(data.shape[1])
    transform = np.fft.fftshift(np.fft.fft2(data * wy[:, None] * wx[None, :]))
    qx = np.fft.fftshift(np.fft.fftfreq(data.shape[1], d=dx_nm))
    qy = np.fft.fftshift(np.fft.fftfreq(data.shape[0], d=dy_nm))
    return qx, qy, np.abs(transform)


def refine_peak(
    image: np.ndarray,
    dx_nm: float,
    dy_nm: float,
    seed: tuple[float, float],
    radius: float = 0.11,
) -> dict:
    qx, qy, amplitude = windowed_fft(image, dx_nm, dy_nm)
    qx_grid, qy_grid = np.meshgrid(qx, qy)
    region = np.hypot(qx_grid - seed[0], qy_grid - seed[1]) <= radius
    local = amplitude[region]
    floor = np.nanpercentile(local, 25)
    weights = np.maximum(local - floor, 0.0) ** 2
    if not np.any(weights > 0):
        raise RuntimeError(f"No Fourier weight near seed {seed}")
    x_value = float(np.sum(weights * qx_grid[region]) / np.sum(weights))
    y_value = float(np.sum(weights * qy_grid[region]) / np.sum(weights))
    return {
        "qx_cycles_per_nm": x_value,
        "qy_cycles_per_nm": y_value,
        "magnitude_cycles_per_nm": float(np.hypot(x_value, y_value)),
        "fft_resolution_cycles_per_nm": [float(abs(qx[1] - qx[0])), float(abs(qy[1] - qy[0]))],
    }


def demodulated_amplitude(
    image: np.ndarray,
    qx_cycles_per_nm: float,
    qy_cycles_per_nm: float,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    sigma_nm: float,
) -> np.ndarray:
    phase = np.exp(
        -2j
        * np.pi
        * (qy_cycles_per_nm * y_nm[:, None] + qx_cycles_per_nm * x_nm[None, :])
    )
    demodulated = image * phase
    dx = float(np.median(np.diff(x_nm)))
    dy = float(np.median(np.diff(y_nm)))
    sigma_pixels = (sigma_nm / dy, sigma_nm / dx)
    real = ndimage.gaussian_filter(demodulated.real, sigma=sigma_pixels, mode="reflect")
    imag = ndimage.gaussian_filter(demodulated.imag, sigma=sigma_pixels, mode="reflect")
    return np.hypot(real, imag)


def derive_boundary(
    image: np.ndarray,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    left_vectors: list[dict],
    right_vectors: list[dict],
    amplitude_sigma_nm: float = 1.7,
    derivative_sigma_nm: float = 0.7,
) -> dict:
    left_components = [
        demodulated_amplitude(
            image,
            item["qx_cycles_per_nm"],
            item["qy_cycles_per_nm"],
            x_nm,
            y_nm,
            amplitude_sigma_nm,
        )
        for item in left_vectors
    ]
    right_components = [
        demodulated_amplitude(
            image,
            item["qx_cycles_per_nm"],
            item["qy_cycles_per_nm"],
            x_nm,
            y_nm,
            amplitude_sigma_nm,
        )
        for item in right_vectors
    ]
    left_amplitude = sum(a / np.nanmedian(a) for a in left_components) / len(left_components)
    right_amplitude = sum(a / np.nanmedian(a) for a in right_components) / len(right_components)
    score = np.log((right_amplitude + 1e-12) / (left_amplitude + 1e-12))
    dx = float(np.median(np.diff(x_nm)))
    dy = float(np.median(np.diff(y_nm)))
    smooth = ndimage.gaussian_filter(score, sigma=(1.2 / dy, 1.2 / dx), mode="reflect")
    derivative = ndimage.gaussian_filter1d(
        smooth, sigma=derivative_sigma_nm / dx, axis=1, order=1, mode="reflect"
    )
    x_search = (x_nm >= 12.0) & (x_nm <= 37.0)
    candidates = np.where(x_search)[0]
    row_x = np.full(len(y_nm), np.nan)
    row_strength = np.full(len(y_nm), np.nan)
    for row in range(len(y_nm)):
        index = candidates[np.argmax(derivative[row, candidates])]
        row_x[row] = x_nm[index]
        row_strength[row] = derivative[row, index]
    valid = (y_nm >= 2.5) & (y_nm <= 37.5)
    valid &= row_strength > np.nanpercentile(row_strength[valid], 20)
    coefficients = np.polyfit(y_nm[valid], row_x[valid], 2)
    for _ in range(3):
        residual = row_x - np.polyval(coefficients, y_nm)
        center = np.nanmedian(residual[valid])
        scale = 1.4826 * np.nanmedian(np.abs(residual[valid] - center))
        valid &= np.abs(residual - center) <= max(2.5 * scale, 1.0)
        coefficients = np.polyfit(y_nm[valid], row_x[valid], 2)
    boundary_x = np.polyval(coefficients, y_nm)
    residual = row_x[valid] - boundary_x[valid]
    uncertainty_nm = float(
        1.4826 * np.nanmedian(np.abs(residual - np.nanmedian(residual)))
    )
    distance_nm = x_nm[None, :] - boundary_x[:, None]
    return {
        "left_amplitude": left_amplitude,
        "right_amplitude": right_amplitude,
        "score": score,
        "smooth_score": smooth,
        "boundary_x_nm": boundary_x,
        "distance_nm": distance_nm,
        "coefficients": coefficients,
        "uncertainty_nm": uncertainty_nm,
        "valid_rows": valid,
        "row_x_nm": row_x,
    }


def local_rms(
    image: np.ndarray,
    dx_nm: float,
    dy_nm: float,
    highpass_nm: float = 1.8,
    rms_nm: float = 0.7,
) -> np.ndarray:
    flat = flatten_linewise(image)
    low = ndimage.gaussian_filter(
        flat, sigma=(highpass_nm / dy_nm, highpass_nm / dx_nm), mode="reflect"
    )
    residual = flat - low
    return np.sqrt(
        ndimage.gaussian_filter(
            residual * residual, sigma=(rms_nm / dy_nm, rms_nm / dx_nm), mode="reflect"
        )
    )


def profile_and_fwhm(local_amplitude: np.ndarray, distance_nm: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    bins = np.linspace(-15.0, 15.0, 61)
    centers = 0.5 * (bins[:-1] + bins[1:])
    profile = np.array(
        [
            np.nanmedian(local_amplitude[(distance_nm >= lower) & (distance_nm < upper)])
            for lower, upper in zip(bins[:-1], bins[1:])
        ]
    )
    smooth = ndimage.gaussian_filter1d(profile, sigma=1.0, mode="nearest")
    baseline = float(np.nanmedian(smooth[np.abs(centers) >= 8.0]))
    central = np.abs(centers) <= 3.0
    peak_index = np.where(central)[0][np.argmax(smooth[central])]
    peak = float(smooth[peak_index])
    half_level = baseline + 0.5 * (peak - baseline)
    above = np.where(smooth >= half_level)[0]
    groups = np.split(above, np.where(np.diff(above) > 1)[0] + 1)
    group = next((item for item in groups if peak_index in item), np.array([peak_index]))
    left = float(centers[group[0]] - 0.25)
    right = float(centers[group[-1]] + 0.25)
    return centers, profile, max(right - left, 0.5)


def analyze_bias_map(
    image: np.ndarray,
    distance_nm: np.ndarray,
    dx_nm: float,
    dy_nm: float,
    highpass_nm: float = 1.8,
    rms_nm: float = 0.7,
) -> dict:
    amplitude = local_rms(image, dx_nm, dy_nm, highpass_nm, rms_nm)
    boundary_mask = np.abs(distance_nm) <= 3.0
    bulk_mask = (np.abs(distance_nm) >= 8.0) & (np.abs(distance_nm) <= 15.0)
    boundary_value = float(np.nanmedian(amplitude[boundary_mask]))
    bulk_values = amplitude[bulk_mask]
    bulk_value = float(np.nanmedian(bulk_values))
    bulk_mad = float(1.4826 * np.nanmedian(np.abs(bulk_values - bulk_value)))
    anomaly_z = (boundary_value - bulk_value) / max(bulk_mad, 1e-30)
    centers, profile, width_nm = profile_and_fwhm(amplitude, distance_nm)
    return {
        "local_amplitude": amplitude,
        "profile_distance_nm": centers,
        "profile": profile,
        "boundary_median_A": boundary_value,
        "bulk_median_A": bulk_value,
        "boundary_excess_z": float(anomaly_z),
        "spatial_fwhm_nm": float(width_nm),
        "detected": bool(anomaly_z > 1.0),
    }


def straighten(
    image: np.ndarray,
    boundary_x_nm: np.ndarray,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    half_width_nm: float = 8.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dx = float(np.median(np.diff(x_nm)))
    u_nm = np.arange(-half_width_nm, half_width_nm + 0.5 * dx, dx)
    strip = np.empty((len(y_nm), len(u_nm)), dtype=float)
    for row in range(len(y_nm)):
        strip[row] = np.interp(
            boundary_x_nm[row] + u_nm, x_nm, image[row], left=np.nan, right=np.nan
        )
    valid = np.isfinite(strip).all(axis=1) & (y_nm >= 2.5) & (y_nm <= 37.5)
    return u_nm, y_nm[valid], strip[valid]


def strip_spectrum(
    strip: np.ndarray, du_nm: float, dy_nm: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = strip.copy()
    data -= np.nanmean(data, axis=1, keepdims=True)
    data -= ndimage.gaussian_filter(data, sigma=(20, 20), mode="reflect")
    wy = np.hanning(data.shape[0])
    wx = np.hanning(data.shape[1])
    transform = np.fft.fftshift(np.fft.fft2(data * wy[:, None] * wx[None, :]))
    q_u = np.fft.fftshift(np.fft.fftfreq(data.shape[1], d=du_nm))
    q_y = np.fft.fftshift(np.fft.fftfreq(data.shape[0], d=dy_nm))
    return q_u, q_y, np.abs(transform)


def analyze_fourier_excess(
    map_images: list[np.ndarray],
    active_indices: list[int],
    control_indices: list[int],
    boundary_x_nm: np.ndarray,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    half_width_nm: float = 8.0,
) -> dict:
    spectra = []
    strips = []
    for image in map_images:
        u_nm, strip_y_nm, strip = straighten(
            flatten_linewise(image), boundary_x_nm, x_nm, y_nm, half_width_nm
        )
        q_u, q_y, amplitude = strip_spectrum(
            strip, float(np.median(np.diff(u_nm))), float(np.median(np.diff(strip_y_nm)))
        )
        strips.append(strip)
        spectra.append(amplitude)
    q_u_grid, q_y_grid = np.meshgrid(q_u, q_y)
    radius = np.hypot(q_u_grid, q_y_grid)
    annulus = (radius >= 0.15) & (radius <= 1.25)
    normalized = [amplitude / np.nanmedian(amplitude[annulus]) for amplitude in spectra]
    active = np.mean([normalized[index] for index in active_indices], axis=0)
    control = np.mean([normalized[index] for index in control_indices], axis=0)
    log_ratio = ndimage.gaussian_filter(
        np.log((active + 1e-12) / (control + 1e-12)), sigma=1.0
    )
    search = annulus & (q_y_grid > 0.05) & (np.abs(q_u_grid) > 0.05)
    peak_flat = int(np.argmax(np.where(search, log_ratio, -np.inf)))
    iy, ix = np.unravel_index(peak_flat, log_ratio.shape)
    peak_u = float(q_u[ix])
    peak_y = float(q_y[iy])
    neighborhood = np.hypot(q_u_grid - peak_u, q_y_grid - peak_y) <= 0.20
    local_floor = float(np.nanpercentile(log_ratio[neighborhood], 20))
    weights = np.maximum(log_ratio[neighborhood] - local_floor, 0.0)
    coordinates = np.column_stack([q_u_grid[neighborhood], q_y_grid[neighborhood]])
    center = np.average(coordinates, axis=0, weights=weights)
    centered = coordinates - center
    covariance = (centered * weights[:, None]).T @ centered / max(np.sum(weights), 1e-30)
    eigenvalues = np.sort(np.maximum(np.linalg.eigvalsh(covariance), 0.0))[::-1]
    fwhm_axes = 2.355 * np.sqrt(eigenvalues)
    return {
        "q_u": q_u,
        "q_y": q_y,
        "log_ratio": log_ratio,
        "strips": strips,
        "feature_q_u_cycles_per_nm": peak_u,
        "feature_q_y_cycles_per_nm": peak_y,
        "feature_q_magnitude_cycles_per_nm": float(np.hypot(peak_u, peak_y)),
        "feature_angle_deg_from_u": float(np.degrees(np.arctan2(peak_y, peak_u))),
        "feature_log_power_ratio": float(log_ratio[iy, ix]),
        "feature_fwhm_major_cycles_per_nm": float(fwhm_axes[0]),
        "feature_fwhm_minor_cycles_per_nm": float(fwhm_axes[1]),
    }


def read_spectra(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    values = np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8")
    bias = values["bias_mV"].astype(float)
    curves = {
        point: values[f"{point}_didv_A"].astype(float)
        for point in ("P01", "P02", "P03")
    }
    return bias, curves


def spectral_candidate(
    bias_mV: np.ndarray,
    curve: np.ndarray,
    references: list[np.ndarray],
    window: int = 21,
) -> dict:
    delta = (curve - np.mean(references, axis=0)) * 1e12
    smooth = savgol_filter(delta, window, 3)
    baseline_mask = ((bias_mV >= 0) & (bias_mV <= 80)) | (
        (bias_mV >= 380) & (bias_mV <= 500)
    )
    baseline_coefficients = np.polyfit(bias_mV[baseline_mask], smooth[baseline_mask], 1)
    baseline = np.polyval(baseline_coefficients, bias_mV)
    corrected = smooth - baseline
    search = (bias_mV >= 80) & (bias_mV <= 360)
    peak_index = np.where(search)[0][np.argmax(corrected[search])]
    peak_height = float(corrected[peak_index])
    half = peak_height / 2.0
    above = np.where((corrected >= half) & search)[0]
    groups = np.split(above, np.where(np.diff(above) > 1)[0] + 1)
    group = next((item for item in groups if peak_index in item), np.array([peak_index]))
    left_index, right_index = int(group[0]), int(group[-1])
    return {
        "delta_pA": delta,
        "smooth_pA": smooth,
        "baseline_pA": baseline,
        "corrected_pA": corrected,
        "peak_index": int(peak_index),
        "feature_peak_mV": float(bias_mV[peak_index]),
        "feature_left_edge_mV": float(bias_mV[left_index]),
        "feature_right_edge_mV": float(bias_mV[right_index]),
        "feature_width_mV": float(bias_mV[right_index] - bias_mV[left_index]),
        "feature_midpoint_mV": float(0.5 * (bias_mV[right_index] + bias_mV[left_index])),
        "feature_peak_height_pA": peak_height,
        "smoothing_window_points": int(window),
    }


def analyze_spectra(path: Path, window: int = 21) -> dict:
    bias, curves = read_spectra(path)
    candidates = {}
    point_ids = sorted(curves)
    for point in point_ids:
        references = [curves[other] for other in point_ids if other != point]
        candidates[point] = spectral_candidate(bias, curves[point], references, window)
    selected = max(point_ids, key=lambda point: candidates[point]["feature_peak_height_pA"])
    return {"bias_mV": bias, "curves_A": curves, "candidates": candidates, "selected_point": selected}


def json_ready(value):
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIR)
    parser.add_argument("--skip-figures", action="store_true")
    args = parser.parse_args()
    start = time.perf_counter()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = json.loads((args.data_dir / "metadata.json").read_text(encoding="utf-8"))

    topography = load_frame(args.data_dir / metadata["topography"]["file"])
    x_nm = topography["x_nm"].astype(float)
    y_nm = topography["y_nm"].astype(float)
    dx_nm = float(np.median(np.diff(x_nm)))
    dy_nm = float(np.median(np.diff(y_nm)))
    topo_mean = flatten_linewise(topography["mean"])

    left_slice = slice(0, 220)
    right_slice = slice(300, 512)
    common_vector = refine_peak(topo_mean, dx_nm, dy_nm, (0.650, 0.275), radius=0.10)
    left_vectors = [
        refine_peak(topo_mean[:, left_slice], dx_nm, dy_nm, (0.698, -0.825), radius=0.12),
        refine_peak(topo_mean[:, left_slice], dx_nm, dy_nm, (0.233, 1.000), radius=0.12),
    ]
    right_vectors = [
        refine_peak(topo_mean[:, right_slice], dx_nm, dy_nm, (0.121, 0.500), radius=0.12),
        refine_peak(topo_mean[:, right_slice], dx_nm, dy_nm, (0.423, -0.425), radius=0.12),
        refine_peak(topo_mean[:, right_slice], dx_nm, dy_nm, (0.302, 0.975), radius=0.12),
    ]
    boundary = derive_boundary(topo_mean, x_nm, y_nm, left_vectors, right_vectors)
    boundary_sample_y = np.array([5.0, 20.0, 35.0])
    boundary_sample_x = np.polyval(boundary["coefficients"], boundary_sample_y)

    map_frames = []
    map_results = []
    for item in metadata["maps"]:
        frame = load_frame(args.data_dir / item["file"])
        map_frames.append(frame)
        analysis = analyze_bias_map(frame["mean"], boundary["distance_nm"], dx_nm, dy_nm)
        map_results.append(
            {
                "id": item["id"],
                "bias_mV": float(item["bias_mV"]),
                "boundary_median_pA": analysis["boundary_median_A"] * 1e12,
                "bulk_median_pA": analysis["bulk_median_A"] * 1e12,
                "boundary_excess_z": analysis["boundary_excess_z"],
                "spatial_fwhm_nm": analysis["spatial_fwhm_nm"],
                "detected": analysis["detected"],
                "registration_shift_pixels_yx": frame["registration_shift_pixels_yx"],
                "forward_retrace_correlation": frame["forward_retrace_correlation"],
                "_analysis": analysis,
            }
        )
    active_indices = [index for index, item in enumerate(map_results) if item["detected"]]
    control_indices = [index for index, item in enumerate(map_results) if not item["detected"]]
    fourier = analyze_fourier_excess(
        [frame["mean"] for frame in map_frames],
        active_indices,
        control_indices,
        boundary["boundary_x_nm"],
        x_nm,
        y_nm,
    )
    spectra = analyze_spectra(args.data_dir / metadata["spectra"]["file"], window=21)
    selected_spectrum = spectra["candidates"][spectra["selected_point"]]

    forward_boundary = derive_boundary(
        flatten_linewise(topography["forward"]), x_nm, y_nm, left_vectors, right_vectors
    )
    retrace_boundary = derive_boundary(
        flatten_linewise(topography["registered_retrace"]), x_nm, y_nm, left_vectors, right_vectors
    )
    validation_boundary_rmse = float(
        np.sqrt(np.mean((forward_boundary["boundary_x_nm"] - retrace_boundary["boundary_x_nm"]) ** 2))
    )
    forward_detected = []
    retrace_detected = []
    validation_rows = []
    for index, (item, frame) in enumerate(zip(metadata["maps"], map_frames)):
        forward_result = analyze_bias_map(
            frame["forward"], boundary["distance_nm"], dx_nm, dy_nm
        )
        retrace_result = analyze_bias_map(
            frame["registered_retrace"], boundary["distance_nm"], dx_nm, dy_nm
        )
        if forward_result["detected"]:
            forward_detected.append(float(item["bias_mV"]))
        if retrace_result["detected"]:
            retrace_detected.append(float(item["bias_mV"]))
        validation_rows.append(
            {
                "id": item["id"],
                "bias_mV": float(item["bias_mV"]),
                "forward_excess_z": forward_result["boundary_excess_z"],
                "retrace_excess_z": retrace_result["boundary_excess_z"],
            }
        )

    robustness_boundaries = []
    robustness_bias = []
    for amplitude_sigma, derivative_sigma in [(1.4, 0.5), (1.7, 0.7), (2.0, 0.9)]:
        trial_boundary = derive_boundary(
            topo_mean,
            x_nm,
            y_nm,
            left_vectors,
            right_vectors,
            amplitude_sigma_nm=amplitude_sigma,
            derivative_sigma_nm=derivative_sigma,
        )
        sample_x = np.polyval(trial_boundary["coefficients"], boundary_sample_y)
        robustness_boundaries.append(sample_x)
        for highpass, rms in [(1.5, 0.6), (1.8, 0.7), (2.1, 0.8)]:
            detected_biases = []
            trial_z = []
            trial_widths = []
            for item, frame in zip(metadata["maps"], map_frames):
                trial = analyze_bias_map(
                    frame["mean"], trial_boundary["distance_nm"], dx_nm, dy_nm, highpass, rms
                )
                trial_z.append(trial["boundary_excess_z"])
                trial_widths.append(trial["spatial_fwhm_nm"])
                if trial["detected"]:
                    detected_biases.append(float(item["bias_mV"]))
            robustness_bias.append(
                {
                    "amplitude_sigma_nm": amplitude_sigma,
                    "derivative_sigma_nm": derivative_sigma,
                    "highpass_nm": highpass,
                    "rms_nm": rms,
                    "detected_biases_mV": detected_biases,
                    "z_scores": trial_z,
                    "widths_nm": trial_widths,
                }
            )
    robustness_boundaries_array = np.asarray(robustness_boundaries)
    fourier_trials = [
        analyze_fourier_excess(
            [frame["mean"] for frame in map_frames],
            active_indices,
            control_indices,
            boundary["boundary_x_nm"],
            x_nm,
            y_nm,
            half_width_nm=half_width,
        )
        for half_width in (7.0, 8.0, 9.0)
    ]
    spectral_trials = [
        analyze_spectra(args.data_dir / metadata["spectra"]["file"], window=window)
        for window in (17, 21, 25, 31)
    ]
    spectral_selected = [trial["candidates"][trial["selected_point"]] for trial in spectral_trials]
    fourier_q_trials = np.array(
        [
            [trial["feature_q_u_cycles_per_nm"], trial["feature_q_y_cycles_per_nm"]]
            for trial in fourier_trials
        ]
    )
    boundary_spread = np.ptp(robustness_boundaries_array, axis=0)
    fourier_spread = np.ptp(fourier_q_trials, axis=0)
    spectral_peak_spread = np.ptp([trial["feature_peak_mV"] for trial in spectral_selected])
    spectral_edge_spread = max(
        np.ptp([trial["feature_left_edge_mV"] for trial in spectral_selected]),
        np.ptp([trial["feature_right_edge_mV"] for trial in spectral_selected]),
    )
    active_biases = [float(map_results[index]["bias_mV"]) for index in active_indices]
    active_widths = [map_results[index]["spatial_fwhm_nm"] for index in active_indices]
    width_trials = [
        row["widths_nm"][index]
        for row in robustness_bias
        for index in active_indices
    ]
    tolerances = {
        "modulation_vector_component_cycles_per_nm": 0.10,
        "boundary_position_nm": float(max(2.0, np.max(boundary_spread) + 0.5)),
        "active_bias_membership": "exact supplied-bias set; no interpolation between unsampled energies",
        "spatial_fwhm_nm": float(max(3.0, 0.5 * np.ptp(width_trials) + dx_nm)),
        "fourier_feature_component_cycles_per_nm": float(max(0.12, np.max(fourier_spread) + 0.05)),
        "spectral_peak_mV": float(max(15.0, spectral_peak_spread + 5.0)),
        "spectral_edge_mV": float(max(15.0, spectral_edge_spread + 5.0)),
        "spectral_width_mV": 25.0,
    }

    results = {
        "schema_version": "1.0",
        "reciprocal_convention": "cycles/nm; multiply by 2*pi for rad/nm",
        "pixel_size_nm": [dx_nm, dy_nm],
        "modulation_vectors": {
            "common": common_vector,
            "left_domain_selective": left_vectors,
            "right_domain_selective": right_vectors,
        },
        "boundary": {
            "model": "x_nm = a*y_nm^2 + b*y_nm + c",
            "coefficients_a_b_c": boundary["coefficients"].tolist(),
            "sample_y_nm": boundary_sample_y.tolist(),
            "sample_x_nm": boundary_sample_x.tolist(),
            "fit_uncertainty_nm": boundary["uncertainty_nm"],
        },
        "bias_dependence": {
            "decision_threshold_boundary_excess_z": 1.0,
            "active_sampled_biases_mV": active_biases,
            "active_sampled_range_mV": [min(active_biases), max(active_biases)],
            "maps": [
                {key: value for key, value in item.items() if key != "_analysis"}
                for item in map_results
            ],
            "active_spatial_fwhm_nm": active_widths,
        },
        "boundary_fourier_feature": {
            key: value
            for key, value in fourier.items()
            if isinstance(value, (float, int, str, bool))
        },
        "spectroscopy": {
            "selected_point": spectra["selected_point"],
            **{
                key: value
                for key, value in selected_spectrum.items()
                if isinstance(value, (float, int, str, bool))
            },
            "energy_step_mV": float(np.median(np.diff(spectra["bias_mV"]))),
        },
        "independent_retrace_validation": {
            "forward_detected_biases_mV": forward_detected,
            "retrace_detected_biases_mV": retrace_detected,
            "topography_boundary_forward_retrace_rmse_nm": validation_boundary_rmse,
            "maps": validation_rows,
        },
        "tolerances": tolerances,
        "environment": {
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    robustness = {
        "boundary_sample_x_nm_min": np.min(robustness_boundaries_array, axis=0).tolist(),
        "boundary_sample_x_nm_max": np.max(robustness_boundaries_array, axis=0).tolist(),
        "bias_parameter_trials": robustness_bias,
        "fourier_half_width_nm_trials": [7.0, 8.0, 9.0],
        "fourier_q_u_q_y_trials": fourier_q_trials.tolist(),
        "spectral_window_trials": [17, 21, 25, 31],
        "spectral_selected_point_trials": [trial["selected_point"] for trial in spectral_trials],
        "spectral_peak_mV_trials": [trial["feature_peak_mV"] for trial in spectral_selected],
        "spectral_edges_mV_trials": [
            [trial["feature_left_edge_mV"], trial["feature_right_edge_mV"]]
            for trial in spectral_selected
        ],
    }

    if not args.skip_figures:
        plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "figure.dpi": 120})
        qx, qy, fft_amplitude = windowed_fft(topo_mean, dx_nm, dy_nm)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        lo, hi = robust_limits(topo_mean * 1e12)
        axes[0].imshow(
            topo_mean * 1e12,
            origin="lower",
            cmap="gray",
            extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
            vmin=lo,
            vmax=hi,
        )
        axes[0].set_title("Flattened topography")
        q_mask_x = np.abs(qx) <= 1.4
        q_mask_y = np.abs(qy) <= 1.4
        crop = np.log1p(fft_amplitude[np.ix_(q_mask_y, q_mask_x)])
        flo, fhi = robust_limits(crop, (5, 99.8))
        axes[1].imshow(
            crop,
            origin="lower",
            cmap="magma",
            extent=[qx[q_mask_x][0], qx[q_mask_x][-1], qy[q_mask_y][0], qy[q_mask_y][-1]],
            vmin=flo,
            vmax=fhi,
        )
        vector_groups = [(common_vector, "C", "#00FFFF")]
        vector_groups += [(item, f"L{index+1}", "#66C2A5") for index, item in enumerate(left_vectors)]
        vector_groups += [(item, f"R{index+1}", "#FC8D62") for index, item in enumerate(right_vectors)]
        for vector, label, color in vector_groups:
            axes[1].plot(vector["qx_cycles_per_nm"], vector["qy_cycles_per_nm"], "o", color=color, ms=6)
            axes[1].text(vector["qx_cycles_per_nm"] + 0.03, vector["qy_cycles_per_nm"] + 0.03, label, color=color)
        axes[0].set_xlabel("x (nm)")
        axes[0].set_ylabel("y (nm)")
        axes[1].set_xlabel("q_x (cycles/nm)")
        axes[1].set_ylabel("q_y (cycles/nm)")
        axes[1].set_title("Low-q Fourier components")
        fig.savefig(args.output_dir / "reciprocal_space_reference.png", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
        axes[0].imshow(boundary["left_amplitude"], origin="lower", cmap="viridis", extent=[0, 40, 0, 40])
        axes[0].set_title("Left-selective local amplitude")
        axes[1].imshow(boundary["right_amplitude"], origin="lower", cmap="viridis", extent=[0, 40, 0, 40])
        axes[1].set_title("Right-selective local amplitude")
        axes[2].imshow(boundary["smooth_score"], origin="lower", cmap="coolwarm", extent=[0, 40, 0, 40])
        axes[2].set_title("Log amplitude ratio and boundary")
        for ax in axes:
            ax.plot(boundary["boundary_x_nm"], y_nm, color="white", lw=1.8)
            ax.set_xlabel("x (nm)")
            ax.set_ylabel("y (nm)")
        fig.savefig(args.output_dir / "domain_segmentation_reference.png", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(2, 4, figsize=(18, 8), constrained_layout=True)
        for column, item in enumerate(map_results):
            analysis = item["_analysis"]
            lo, hi = robust_limits(analysis["local_amplitude"] * 1e12, (2, 98))
            axes[0, column].imshow(
                analysis["local_amplitude"] * 1e12,
                origin="lower",
                cmap="magma",
                extent=[0, 40, 0, 40],
                vmin=lo,
                vmax=hi,
            )
            axes[0, column].plot(boundary["boundary_x_nm"], y_nm, color="#00FFFF", lw=1.3)
            axes[0, column].set_title(f"{item['bias_mV']:.0f} mV | z={item['boundary_excess_z']:.2f}")
            axes[1, column].plot(
                analysis["profile_distance_nm"], analysis["profile"] * 1e12, color="#2C7BB6"
            )
            axes[1, column].axvspan(-3, 3, color="#FDAE61", alpha=0.25)
            axes[1, column].set_title(f"spatial FWHM={item['spatial_fwhm_nm']:.1f} nm")
            axes[1, column].set_xlabel("signed distance (nm)")
            axes[1, column].set_ylabel("local RMS (pA)")
        fig.savefig(args.output_dir / "bias_dependence_reference.png", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        active_strip = np.mean([fourier["strips"][index] for index in active_indices], axis=0)
        slo, shi = robust_limits(active_strip * 1e12)
        axes[0].imshow(
            active_strip * 1e12,
            origin="lower",
            aspect="auto",
            cmap="RdBu_r",
            extent=[-8.0, 8.0, 2.5, 37.5],
            vmin=slo,
            vmax=shi,
        )
        axes[0].set_title("Mean straightened active-bias strip")
        axes[0].set_xlabel("boundary-normal offset u (nm)")
        axes[0].set_ylabel("analysis y (nm)")
        q_u = fourier["q_u"]
        q_y = fourier["q_y"]
        mask_u = np.abs(q_u) <= 1.3
        mask_y = np.abs(q_y) <= 1.3
        shown = fourier["log_ratio"][np.ix_(mask_y, mask_u)]
        limit = float(np.nanpercentile(np.abs(shown), 98))
        axes[1].imshow(
            shown,
            origin="lower",
            cmap="coolwarm",
            extent=[q_u[mask_u][0], q_u[mask_u][-1], q_y[mask_y][0], q_y[mask_y][-1]],
            vmin=-limit,
            vmax=limit,
        )
        axes[1].plot(
            fourier["feature_q_u_cycles_per_nm"],
            fourier["feature_q_y_cycles_per_nm"],
            marker="x",
            color="yellow",
            ms=10,
            mew=2,
        )
        axes[1].set_xlabel("q_u (cycles/nm)")
        axes[1].set_ylabel("q_y (cycles/nm)")
        axes[1].set_title("Active/control Fourier log-power ratio")
        fig.savefig(args.output_dir / "boundary_feature_fft_reference.png", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True, sharex=True)
        for point, curve in spectra["curves_A"].items():
            axes[0].plot(spectra["bias_mV"], curve * 1e12, label=point)
        axes[0].legend()
        axes[0].set_ylabel("dI/dV proxy (pA)")
        axes[0].set_title("Neutral point spectra")
        axes[1].plot(spectra["bias_mV"], selected_spectrum["delta_pA"], color="0.7", label="raw difference")
        axes[1].plot(spectra["bias_mV"], selected_spectrum["corrected_pA"], color="#008837", label="baseline-corrected")
        axes[1].axvspan(
            selected_spectrum["feature_left_edge_mV"],
            selected_spectrum["feature_right_edge_mV"],
            color="#FDAE61",
            alpha=0.3,
        )
        axes[1].axvline(selected_spectrum["feature_peak_mV"], color="#D73027", ls="--")
        axes[1].legend()
        axes[1].set_xlabel("Bias (mV)")
        axes[1].set_ylabel("difference (pA)")
        axes[1].set_title(
            f"selected {spectra['selected_point']}: peak={selected_spectrum['feature_peak_mV']:.1f} mV, "
            f"FWHM={selected_spectrum['feature_width_mV']:.1f} mV"
        )
        fig.savefig(args.output_dir / "spectroscopy_reference.png", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
        axes[0].plot(forward_boundary["boundary_x_nm"], y_nm, label="forward")
        axes[0].plot(retrace_boundary["boundary_x_nm"], y_nm, label="retrace")
        axes[0].set_xlabel("boundary x (nm)")
        axes[0].set_ylabel("y (nm)")
        axes[0].set_title(f"boundary RMSE={validation_boundary_rmse:.2f} nm")
        axes[0].legend()
        biases = [row["bias_mV"] for row in validation_rows]
        axes[1].plot(biases, [row["forward_excess_z"] for row in validation_rows], "o-", label="forward")
        axes[1].plot(biases, [row["retrace_excess_z"] for row in validation_rows], "s-", label="retrace")
        axes[1].axhline(1.0, color="black", ls="--", label="decision threshold")
        axes[1].set_xlabel("bias (mV)")
        axes[1].set_ylabel("boundary-excess z")
        axes[1].set_title("Independent scan-direction validation")
        axes[1].legend()
        fig.savefig(args.output_dir / "validation_reference.png", dpi=180)
        plt.close(fig)

    elapsed = time.perf_counter() - start
    results["runtime_seconds_this_run"] = float(elapsed)
    (args.output_dir / "reference_results.json").write_text(
        json.dumps(json_ready(results), ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (args.output_dir / "robustness_summary.json").write_text(
        json.dumps(json_ready(robustness), ensure_ascii=False, indent=2), encoding="utf-8"
    )
    with (args.output_dir / "reference_results.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["quantity", "value", "unit_or_note"])
        writer.writerow(["active_sampled_biases", ";".join(f"{v:.0f}" for v in active_biases), "mV"])
        writer.writerow(["boundary_x_at_y_5_20_35", ";".join(f"{v:.4f}" for v in boundary_sample_x), "nm"])
        writer.writerow(["fourier_feature_q_u", fourier["feature_q_u_cycles_per_nm"], "cycles/nm"])
        writer.writerow(["fourier_feature_q_y", fourier["feature_q_y_cycles_per_nm"], "cycles/nm"])
        writer.writerow(["spectral_selected_point", spectra["selected_point"], "neutral point ID"])
        writer.writerow(["spectral_peak", selected_spectrum["feature_peak_mV"], "mV"])
        writer.writerow(["spectral_left_edge", selected_spectrum["feature_left_edge_mV"], "mV"])
        writer.writerow(["spectral_right_edge", selected_spectrum["feature_right_edge_mV"], "mV"])
        writer.writerow(["spectral_width", selected_spectrum["feature_width_mV"], "mV"])
        writer.writerow(["runtime", elapsed, "seconds"])
    print(json.dumps({
        "runtime_seconds": elapsed,
        "active_sampled_biases_mV": active_biases,
        "selected_spectral_point": spectra["selected_point"],
        "output_dir": str(args.output_dir),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
