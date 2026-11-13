#!/usr/bin/env python3
"""Deterministic Harbor verifier for the GdTe2 boundary-charge-order STM/STS task.

Scores the agent's deliverables in the submission directory against the
result-only rubric in grading_no_tool.md (10 checkpoints, 100 points total).
Everything is computed from the submitted summary.json / analysis_report.md /
PNG files; no network and no LLM judge are required.

Reference numbers below are the values recomputed by solution/reference_analysis.py
from /app/data (they match answer/reference_results.json).

Figure-content points are awarded by an automated proxy: the PNG must exist, be a
valid non-blank image of reasonable size. A human grader would inspect the actual
annotations; the deterministic verifier cannot, so existence + validity is used.

Writes:
  <logs>/reward.json          multi-metric reward (primary key "reward" in [0, 1])
  <logs>/reward.txt           plain-text normalized reward in [0, 1]
  <logs>/grading_report.json  human-readable per-checkpoint breakdown
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from pathlib import Path

# --------------------------------------------------------------------------- #
# Reference values (recomputed by solution/reference_analysis.py from raw data).
# --------------------------------------------------------------------------- #
REF_COMMON = (0.655, 0.275, 0.710)
REF_LEFT = [(0.686, -0.828, 1.075), (0.221, 1.007, 1.031)]
REF_RIGHT = [(0.141, 0.490, 0.510), (0.422, -0.433, 0.604), (0.286, 0.981, 1.022)]
VEC_TOL = 0.10

REF_BOUNDARY_Y = [5.0, 20.0, 35.0]
REF_BOUNDARY_X = {5.0: 18.30, 20.0: 24.12, 35.0: 29.55}
BOUNDARY_TOL_FULL = 2.0
BOUNDARY_TOL_PARTIAL = 3.5
REF_UNC_LO, REF_UNC_HI = 0.3, 2.5

REF_ACTIVE = {300.0, 400.0}
REF_CONTROLS = {100.0, 900.0}
REF_SAMPLED = {100.0, 300.0, 400.0, 900.0}
REF_Z = {100.0: -0.23, 300.0: 3.14, 400.0: 1.75, 900.0: 0.13}
Z_THRESHOLD = 1.0
REF_FWHM = {300.0: 8.5, 400.0: 9.5}
FWHM_TOL = 3.0

REF_FEAT_QU, REF_FEAT_QY = 0.621, 0.513
FEAT_COMP_TOL = 0.12
REF_FEAT_MAG = 0.806
FEAT_MAG_TOL = 0.15
REF_FEAT_ANGLE = 39.55
FEAT_ANGLE_TOL = 15.0

REF_SPEC_POINT = "P02"
REF_PEAK, REF_LEFT_EDGE, REF_RIGHT_EDGE, REF_WIDTH = 150.0, 127.5, 202.5, 75.0
SPEC_PEAK_TOL = 15.0
SPEC_WIDTH_TOL = 25.0

REF_RMSE_LO, REF_RMSE_HI = 0.0, 3.68

REQUIRED_FIGURES = [
    "reciprocal_space.png",
    "domain_segmentation.png",
    "bias_dependence.png",
    "boundary_feature_fft.png",
    "spectroscopy_comparison.png",
    "validation_retrace.png",
]
REQUIRED_FILES = ["analysis_report.md", "summary.json", "analysis.py"] + REQUIRED_FIGURES
TOP_FIELDS = [
    "schema_version", "reciprocal_convention", "common_modulation_vector",
    "left_domain_vectors", "right_domain_vectors", "boundary", "boundary_anomaly",
    "fourier_feature", "spectroscopy", "retrace_validation",
]


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def _norm(k):
    return re.sub(r"[^a-z0-9]", "", str(k).lower())


def num(x):
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def vget(d, names):
    """Fetch the first matching key from a dict, normalized-case-insensitively."""
    if not isinstance(d, dict):
        return None
    norm = {_norm(k): v for k, v in d.items()}
    for n in names:
        nk = _norm(n)
        if nk in norm:
            return norm[nk]
    return None


def read_text(path):
    try:
        return Path(path).read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""


def has_any(text, keywords):
    low = text.lower()
    return any(k.lower() in low for k in keywords)


def count_topics(text, topic_kw):
    return sum(1 for kws in topic_kw if has_any(text, kws))


def figure_ok(path, min_bytes=4000):
    p = Path(path)
    if not p.exists() or p.stat().st_size < min_bytes:
        return False
    try:
        from PIL import Image
        import numpy as np
        im = Image.open(p)
        im.load()
        if im.width < 64 or im.height < 64:
            return False
        arr = np.asarray(im.convert("L"), dtype=float)
        return bool(arr.std() > 1.0)
    except Exception:
        return p.stat().st_size >= min_bytes


def parse_vec(item):
    """Return (qx, qy, magnitude) from a submitted vector dict (fuzzy keys)."""
    if not isinstance(item, dict):
        return None
    qx = num(vget(item, ["qx", "q_x", "qx_cycles_per_nm", "qx_cycles", "qx_nm"]))
    qy = num(vget(item, ["qy", "q_y", "qy_cycles_per_nm", "qy_cycles", "qy_nm"]))
    mag = num(vget(item, ["magnitude", "mag", "q", "q_magnitude", "magnitude_cycles_per_nm", "abs", "norm"]))
    if qx is None and qy is None and mag is None:
        return None
    if mag is None and qx is not None and qy is not None:
        mag = math.hypot(qx, qy)
    return (qx, qy, mag)


def vec_list(raw):
    out = []
    if isinstance(raw, list):
        for it in raw:
            v = parse_vec(it)
            if v is not None:
                out.append(v)
    return out


def comp_match_count(sub, ref, tol):
    """Best count (0..2) of components within tol under a global sign flip."""
    if sub is None:
        return 0
    best = 0
    for s in (1.0, -1.0):
        c = 0
        if sub[0] is not None and abs(sub[0] - s * ref[0]) <= tol:
            c += 1
        if sub[1] is not None and abs(sub[1] - s * ref[1]) <= tol:
            c += 1
        best = max(best, c)
    return best


def angle_of(v):
    if v is None or v[0] is None or v[1] is None:
        return None
    return math.degrees(math.atan2(v[1], v[0]))


def undirected_angle_diff(a, b):
    if a is None or b is None:
        return None
    d = abs(a - b) % 180.0
    return min(d, 180.0 - d)


def ref_mag(ref):
    return math.hypot(ref[0], ref[1])


def vec_pair_score(sub, ref, tol, angle_tol=15.0):
    """2.5 if both components match; 1.0 if only magnitude or direction; else 0."""
    if sub is None:
        return 0.0
    if comp_match_count(sub, ref, tol) == 2:
        return 2.5
    mag_ok = sub[2] is not None and abs(sub[2] - ref_mag(ref)) <= tol
    adiff = undirected_angle_diff(angle_of(sub), angle_of(ref))
    dir_ok = adiff is not None and adiff <= angle_tol
    if mag_ok or dir_ok:
        return 1.0
    return 0.0


def parse_bias_key(k):
    m = re.search(r"-?\d+\.?\d*", str(k))
    return float(m.group()) if m else None


def bias_dict(raw):
    """Normalize a {bias: value} dict to {float_bias: float_value}."""
    out = {}
    if isinstance(raw, dict):
        for k, v in raw.items():
            f = parse_bias_key(k)
            if f is not None:
                out[f] = num(v)
    elif isinstance(raw, list):
        # tolerate [{"bias_mV":..,"z":..}] style
        for it in raw:
            if isinstance(it, dict):
                b = num(vget(it, ["bias_mV", "bias", "energy_mV", "energy"]))
                val = num(vget(it, ["z", "boundary_excess_z", "excess_z", "score", "value", "fwhm_nm", "spatial_fwhm_nm"]))
                if b is not None:
                    out[b] = val
    return out


def bias_list(raw):
    """Normalize a list of biases to a set of floats."""
    out = set()
    if isinstance(raw, list):
        for it in raw:
            f = num(it) if not isinstance(it, dict) else num(vget(it, ["bias_mV", "bias", "energy_mV"]))
            if f is None:
                f = parse_bias_key(it)
            if f is not None:
                out.add(round(float(f), 1))
    return out


def lookup_bias(d, target, tol=0.6):
    if not d:
        return None
    best, bestd = None, tol + 1
    for k, v in d.items():
        dd = abs(k - target)
        if dd < bestd:
            best, bestd = v, dd
    return best if bestd <= tol else None


def decimals_ok(x, max_dec):
    if x is None:
        return True
    s = repr(float(x))
    if "e" in s or "E" in s:
        return False
    frac = s.split(".")[1] if "." in s else ""
    return len(frac.rstrip("0")) <= max_dec


# --------------------------------------------------------------------------- #
# Checkpoint 1: output completeness & schema (6)
# --------------------------------------------------------------------------- #
def cp1_completeness(sub, summary, report):
    d = {"files": {}, "summary_parsed": summary is not None}
    n = 0
    for name in REQUIRED_FILES:
        ok = (sub / name).exists()
        if name == "summary.json":
            ok = ok and summary is not None
        d["files"][name] = bool(ok)
        n += 1 if ok else 0
    fields_present = 0
    if isinstance(summary, dict):
        for f in TOP_FIELDS:
            if f in summary:
                fields_present += 1
    d["top_fields_present"] = fields_present
    fields_ok = fields_present == len(TOP_FIELDS)
    d["all_top_fields"] = fields_ok

    has_summary = summary is not None
    has_report = (sub / "analysis_report.md").exists() and len(report.strip()) > 0

    if has_summary:
        if n == len(REQUIRED_FILES) and fields_ok:
            earned = 6.0
        elif n == len(REQUIRED_FILES):
            earned = 5.0
        elif n >= len(REQUIRED_FILES) - 1:
            earned = 4.0
        elif n >= len(REQUIRED_FILES) - 2:
            earned = 3.0
        else:
            earned = 3.0
    else:
        if has_report:
            earned = 2.0 if len(report) > 800 else 1.0
        else:
            earned = 0.0
    d["files_present"] = n
    return earned, 6.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 2: units & coordinate conventions (4)
# --------------------------------------------------------------------------- #
def cp2_units(sub, summary, report):
    d = {}
    conv = ""
    if isinstance(summary, dict):
        conv = str(summary.get("reciprocal_convention", ""))
    conv_cycles = "cycles" in conv.lower()
    report_cycles = has_any(report, ["cycles/nm", "cycles per nm"])
    report_rad = has_any(report, ["rad/nm"]) and has_any(report, ["2π", "2*pi", "2\\pi", "two pi"])
    # 2 pts: convention field cycles/nm AND report explains 2pi conversion
    if conv_cycles and report_rad:
        conv_score = 2.0
    elif conv_cycles or report_cycles:
        conv_score = 1.0
    else:
        conv_score = 0.0
    d["convention_cycles"] = conv_cycles
    d["report_rad_conversion"] = report_rad
    d["conv_score"] = conv_score

    # 1 pt: x/y tied to columns/rows (or pixel scale) and lengths in nm
    coord = (has_any(report, ["列", "column", "像素尺度", "pixel size", "pixel_scale", "像素"])) and \
            (has_any(report, ["行", "row", "nm"]))
    coord_score = 1.0 if coord else 0.0
    d["coord_score"] = coord_score

    # 1 pt: sensible rounding in summary (q<=3, length<=2, bias<=1 decimals)
    dec_ok = True
    if isinstance(summary, dict):
        cmv = parse_vec(summary.get("common_modulation_vector"))
        if cmv:
            for v in cmv[:3]:
                if not decimals_ok(v, 3):
                    dec_ok = False
        bnd = summary.get("boundary") if isinstance(summary.get("boundary"), dict) else {}
        for x in (bnd.get("sample_x_nm") or []):
            if not decimals_ok(num(x), 2):
                dec_ok = False
        for b in (summary.get("boundary_anomaly", {}) or {}).get("active_sampled_biases_mV", []) or []:
            if not decimals_ok(num(b), 1):
                dec_ok = False
    else:
        dec_ok = False
    dec_score = 1.0 if dec_ok else 0.0
    d["decimals_ok"] = dec_ok

    return conv_score + coord_score + dec_score, 4.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 3: common modulation vector (6)
# --------------------------------------------------------------------------- #
def cp3_common(sub, summary):
    d = {}
    cmv = parse_vec(vget(summary, ["common_modulation_vector"]) if isinstance(summary, dict) else None)
    ref = REF_COMMON
    cm = comp_match_count(cmv, ref, VEC_TOL)
    comp_score = cm * 2.0  # 2 pts per component, max 4
    mag = cmv[2] if cmv else None
    mag_score = 1.0 if (mag is not None and abs(mag - ref[2]) <= VEC_TOL) else 0.0
    fig_score = 1.0 if figure_ok(sub / "reciprocal_space.png") else 0.0
    d.update({"components_matched": cm, "comp_score": comp_score,
              "magnitude": mag, "mag_score": mag_score, "fig_score": fig_score})
    return comp_score + mag_score + fig_score, 6.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 4: two-sided selective vectors (10)
# --------------------------------------------------------------------------- #
def cp4_selective(sub, summary):
    d = {}
    left = vec_list(vget(summary, ["left_domain_vectors"]) if isinstance(summary, dict) else None)
    right = vec_list(vget(summary, ["right_domain_vectors"]) if isinstance(summary, dict) else None)

    # Left: two reference vectors, 2.5 each (cap 5).
    left_score = 0.0
    left_detail = []
    for ref in REF_LEFT:
        best = max((vec_pair_score(sv, ref, VEC_TOL) for sv in left), default=0.0)
        left_detail.append(round(best, 2))
        left_score += best
    left_score = min(left_score, 5.0)

    # Right: three acceptable references, only top two count (cap 5).
    right_per_ref = []
    for ref in REF_RIGHT:
        best = max((vec_pair_score(sv, ref, VEC_TOL) for sv in right), default=0.0)
        right_per_ref.append(round(best, 2))
    right_score = min(sum(sorted(right_per_ref, reverse=True)[:2]), 5.0)

    d.update({"left_per_ref": left_detail, "left_score": round(left_score, 2),
              "right_per_ref": right_per_ref, "right_score": round(right_score, 2),
              "n_left": len(left), "n_right": len(right)})
    return left_score + right_score, 10.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 5: data-driven boundary localization (14)
# --------------------------------------------------------------------------- #
def cp5_boundary(sub, summary, report):
    d = {}
    bnd = vget(summary, ["boundary"]) if isinstance(summary, dict) else None
    bnd = bnd if isinstance(bnd, dict) else {}
    ys = [num(v) for v in (bnd.get("sample_y_nm") or [])]
    xs = [num(v) for v in (bnd.get("sample_x_nm") or [])]

    pos_score = 0.0
    pos_detail = {}
    for ty in REF_BOUNDARY_Y:
        rx = REF_BOUNDARY_X[ty]
        # locate submitted x nearest to target y
        sx = None
        if ys and xs and len(ys) == len(xs):
            best, bestd = None, 1e9
            for yy, xx in zip(ys, xs):
                if yy is None or xx is None:
                    continue
                dd = abs(yy - ty)
                if dd < bestd:
                    best, bestd = xx, dd
            if bestd <= 1.0:
                sx = best
        elif xs and len(xs) == 3:
            sx = xs[REF_BOUNDARY_Y.index(ty)]
        if sx is None:
            pos_detail[ty] = {"submitted": None, "pts": 0.0}
            continue
        dx = abs(sx - rx)
        pts = 3.0 if dx <= BOUNDARY_TOL_FULL else (1.0 if dx <= BOUNDARY_TOL_PARTIAL else 0.0)
        pos_detail[ty] = {"submitted": round(sx, 3), "ref": rx, "dx": round(dx, 3), "pts": pts}
        pos_score += pts

    # 2 pts: sample_y correct AND curve increasing (lower-left -> upper-right)
    y_ok = len(ys) == 3 and all(a is not None and abs(a - b) <= 1.0 for a, b in zip(sorted(ys), REF_BOUNDARY_Y))
    xvals = [x for x in xs if x is not None]
    monotonic = len(xvals) >= 2 and all(xvals[i] < xvals[i + 1] for i in range(len(xvals) - 1))
    curve_score = 2.0 if (y_ok and monotonic) else (1.0 if (y_ok or monotonic) else 0.0)

    # 2 pts: uncertainty reported within 0.3-2.5 nm
    unc = num(bnd.get("uncertainty_nm"))
    if unc is None:
        unc = num(vget(bnd, ["uncertainty_nm", "uncertainty", "fit_uncertainty_nm", "std_nm"]))
    unc_score = 2.0 if (unc is not None and REF_UNC_LO <= unc <= REF_UNC_HI) else (1.0 if unc is not None else 0.0)

    fig_score = 1.0 if figure_ok(sub / "domain_segmentation.png") else 0.0

    d.update({"positions": pos_detail, "pos_score": round(pos_score, 2),
              "y_ok": y_ok, "monotonic": monotonic, "curve_score": curve_score,
              "uncertainty_nm": unc, "unc_score": unc_score, "fig_score": fig_score})
    return pos_score + curve_score + unc_score + fig_score, 14.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 6: bias-dependent boundary anomaly (16)
# --------------------------------------------------------------------------- #
def cp6_bias_anomaly(sub, summary, report):
    d = {}
    ba = vget(summary, ["boundary_anomaly"]) if isinstance(summary, dict) else None
    ba = ba if isinstance(ba, dict) else {}
    active = bias_list(vget(ba, ["active_sampled_biases_mV", "active_biases_mV", "active_biases"]))

    # Active-set tier (max 10)
    has_both = REF_ACTIVE.issubset(active)
    has_control = bool(active & REF_CONTROLS)
    if active == REF_ACTIVE:
        active_score = 10.0
        tier = "exact"
    elif has_both and has_control:
        active_score = 3.0
        tier = "both_plus_false_control"
    elif active and active.issubset(REF_ACTIVE):
        active_score = 6.0
        tier = "one_correct_only"
    else:
        active_score = 0.0
        tier = "miss"
    d.update({"active": sorted(active), "tier": tier, "active_score": active_score})

    # z-scores (max 4): 300>1, 400>1, 100<1, 900<1
    zraw = bias_dict(vget(ba, ["per_bias_boundary_excess_z", "per_bias_z", "boundary_excess_z"]))
    z_score = 0.0
    z_detail = {}
    for b in sorted(REF_SAMPLED):
        z = lookup_bias(zraw, b)
        if z is None:
            z_detail[b] = {"z": None, "pts": 0.0}
            continue
        expect_active = b in REF_ACTIVE
        ok = (z > Z_THRESHOLD) if expect_active else (z < Z_THRESHOLD)
        pts = 1.0 if ok else 0.0
        z_detail[b] = {"z": round(z, 3), "expected_active": expect_active, "pts": pts}
        z_score += pts
    d["z_detail"] = z_detail
    d["z_score"] = round(z_score, 2)

    # detected consistency (1) - only if the field is actually present, so an
    # empty/partial submission cannot earn it by defaulting to False.
    det = ba.get("detected")
    det_score = 1.0 if ("detected" in ba and bool(det) == (len(active) > 0)) else 0.0
    d["detected"] = det
    d["det_score"] = det_score

    # no fabricated continuous endpoints (1)
    all_sampled = active.issubset(REF_SAMPLED)
    kw = has_any(report, ["已采样", "sampled", "不插值", "不外推", "no interpolation",
                          "no extrapolation", "离散", "discrete", "连续能量"])
    nofab_score = 1.0 if (all_sampled and kw) else 0.0
    d["nofab_score"] = nofab_score

    return active_score + z_score + det_score + nofab_score, 16.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 7: spatial width & boundary Fourier feature (14)
# --------------------------------------------------------------------------- #
def cp7_width_fourier(sub, summary):
    d = {}
    ba = vget(summary, ["boundary_anomaly"]) if isinstance(summary, dict) else None
    ba = ba if isinstance(ba, dict) else {}
    fw = bias_dict(vget(ba, ["active_spatial_fwhm_nm", "spatial_fwhm_nm", "fwhm_nm"]))
    fw_score = 0.0
    fw_detail = {}
    for b in sorted(REF_FWHM):
        v = lookup_bias(fw, b)
        pts = 2.0 if (v is not None and abs(v - REF_FWHM[b]) <= FWHM_TOL) else 0.0
        fw_detail[b] = {"fwhm": v, "ref": REF_FWHM[b], "pts": pts}
        fw_score += pts
    d["fwhm_detail"] = fw_detail
    d["fwhm_score"] = round(fw_score, 2)

    feat = vget(summary, ["fourier_feature"]) if isinstance(summary, dict) else None
    feat = feat if isinstance(feat, dict) else {}
    qu = num(vget(feat, ["q_u_cycles_per_nm", "q_u", "qu", "q_u_cycles"]))
    qy = num(vget(feat, ["q_y_cycles_per_nm", "q_y", "qy", "q_v"]))
    mag = num(vget(feat, ["magnitude_cycles_per_nm", "magnitude", "mag", "q_magnitude"]))
    ang = num(vget(feat, ["angle_deg_from_u", "angle_deg", "angle", "theta_deg"]))
    sub_v = (qu, qy, mag if mag is not None else (math.hypot(qu, qy) if qu is not None and qy is not None else None))
    ref_v = (REF_FEAT_QU, REF_FEAT_QY, REF_FEAT_MAG)
    cm = comp_match_count(sub_v, ref_v, FEAT_COMP_TOL)
    comp_score = cm * 2.0  # q_u, q_y each 2 pts, max 4
    mag_score = 2.0 if (sub_v[2] is not None and abs(sub_v[2] - REF_FEAT_MAG) <= FEAT_MAG_TOL) else 0.0
    adiff = undirected_angle_diff(ang, REF_FEAT_ANGLE)
    if adiff is None:
        adiff = undirected_angle_diff(angle_of(sub_v), REF_FEAT_ANGLE)
    ang_score = 2.0 if (adiff is not None and adiff <= FEAT_ANGLE_TOL) else 0.0
    fig_score = 2.0 if figure_ok(sub / "boundary_feature_fft.png") else 0.0
    d.update({"q_u": qu, "q_y": qy, "comp_score": comp_score, "magnitude": mag,
              "mag_score": mag_score, "angle": ang, "angle_diff": adiff,
              "ang_score": ang_score, "fig_score": fig_score})
    return comp_score + mag_score + ang_score + fig_score + fw_score, 14.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 8: point-spectroscopy local feature (16)
# --------------------------------------------------------------------------- #
def cp8_spectroscopy(sub, summary):
    d = {}
    sp = vget(summary, ["spectroscopy"]) if isinstance(summary, dict) else None
    sp = sp if isinstance(sp, dict) else {}
    pid = vget(sp, ["anomalous_point_id", "selected_point", "point_id", "anomalous_point"])
    pid_score = 4.0 if (pid is not None and str(pid).strip().upper() == REF_SPEC_POINT) else 0.0

    peak = num(vget(sp, ["feature_peak_mV", "peak_mV", "peak", "feature_peak"]))
    if peak is not None and abs(peak - REF_PEAK) <= SPEC_PEAK_TOL:
        peak_score = 3.0
    elif peak is not None and 150.0 <= peak <= 180.0:
        peak_score = 2.0  # equivalent center definition
    else:
        peak_score = 0.0

    left = num(vget(sp, ["left_edge_mV", "feature_left_edge_mV", "left_edge"]))
    right = num(vget(sp, ["right_edge_mV", "feature_right_edge_mV", "right_edge"]))
    left_score = 2.0 if (left is not None and abs(left - REF_LEFT_EDGE) <= SPEC_PEAK_TOL) else 0.0
    right_score = 2.0 if (right is not None and abs(right - REF_RIGHT_EDGE) <= SPEC_PEAK_TOL) else 0.0

    width = num(vget(sp, ["width_mV", "feature_width_mV", "fwhm_mV", "width"]))
    if width is None and left is not None and right is not None:
        width = abs(right - left)
    width_score = 3.0 if (width is not None and abs(width - REF_WIDTH) <= SPEC_WIDTH_TOL) else 0.0

    fig_score = 2.0 if figure_ok(sub / "spectroscopy_comparison.png") else 0.0
    d.update({"point_id": pid, "pid_score": pid_score, "peak": peak, "peak_score": peak_score,
              "left": left, "left_score": left_score, "right": right, "right_score": right_score,
              "width": width, "width_score": width_score, "fig_score": fig_score})
    return pid_score + peak_score + left_score + right_score + width_score + fig_score, 16.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 9: forward/retrace validation (8)
# --------------------------------------------------------------------------- #
def cp9_retrace(sub, summary):
    d = {}
    rv = vget(summary, ["retrace_validation"]) if isinstance(summary, dict) else None
    rv = rv if isinstance(rv, dict) else {}
    fwd = bias_list(vget(rv, ["forward_active_biases_mV", "forward_active", "forward_detected_biases_mV"]))
    ret = bias_list(vget(rv, ["retrace_active_biases_mV", "retrace_active", "retrace_detected_biases_mV"]))

    def dir_score(s):
        if s == REF_ACTIVE:
            return 2.0
        if s and s.issubset(REF_ACTIVE):
            return 1.0
        return 0.0

    fwd_score = dir_score(fwd)
    ret_score = dir_score(ret)

    rmse = num(vget(rv, ["boundary_forward_retrace_rmse_nm", "rmse_nm", "boundary_rmse_nm", "rmse"]))
    rmse_score = 2.0 if (rmse is not None and REF_RMSE_LO <= rmse <= REF_RMSE_HI) else 0.0
    fig_score = 2.0 if figure_ok(sub / "validation_retrace.png") else 0.0
    d.update({"forward": sorted(fwd), "fwd_score": fwd_score, "retrace": sorted(ret),
              "ret_score": ret_score, "rmse": rmse, "rmse_score": rmse_score, "fig_score": fig_score})
    return fwd_score + ret_score + rmse_score + fig_score, 8.0, d


# --------------------------------------------------------------------------- #
# Checkpoint 10: scientific report & tolerance explanation (6)
# --------------------------------------------------------------------------- #
def cp10_report(sub, report):
    d = {}
    topics = [
        ["预处理", "preprocess", "背景", "flatten", "配准", "registration"],
        ["cycles/nm", "fft", "fourier", "傅里叶", "倒空间"],
        ["分区", "partition", "η", "eta", "局域振幅", "segmentation"],
        ["异常分数", "anomaly", "s(e)", "阈值", "threshold", "z 分数", "z-score", "excess"],
        ["fwhm", "宽度", "半高", "width"],
        ["边缘", "edge", "基线", "baseline"],
    ]
    ntopics = count_topics(report, topics)
    def_score = 2.0 if ntopics >= 4 else (1.0 if ntopics >= 2 else 0.0)

    resolution = has_any(report, ["分辨率", "像素尺度", "resolution", "0.078", "间隔", "energy step", "能量步长"])
    uncertainty = has_any(report, ["不确定", "uncertainty", "误差", "鲁棒", "robust", "error"])
    res_score = 2.0 if (resolution and uncertainty) else (1.0 if (resolution or uncertainty) else 0.0)

    paper = has_any(report, ["论文", "paper", "supporting information", "si "])
    rawprio = has_any(report, ["sanity", "重算", "raw", "数据", "ground truth", "不为贴合"])
    raw_score = 1.0 if (paper and rawprio) else 0.0

    retrace = has_any(report, ["retrace", "扫描方向", "scan direction", "反向"])
    secondb = has_any(report, ["第二", "second", "独立", "independent", "扫描方向"])
    scope_score = 1.0 if (retrace and secondb) else 0.0

    d.update({"topics_matched": ntopics, "def_score": def_score,
              "resolution": resolution, "uncertainty": uncertainty, "res_score": res_score,
              "raw_vs_paper": raw_score, "scope": scope_score, "report_len": len(report)})
    return def_score + res_score + raw_score + scope_score, 6.0, d


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def grade(submission_dir: Path, logs_dir: Path):
    sub = Path(submission_dir)
    logs = Path(logs_dir)
    logs.mkdir(parents=True, exist_ok=True)

    summary = None
    try:
        summary = json.loads((sub / "summary.json").read_text(encoding="utf-8"))
    except Exception:
        summary = None
    if not isinstance(summary, dict):
        summary = None
    report = read_text(sub / "analysis_report.md")

    checkpoints = {}
    for key, fn in [
        ("cp1_completeness", lambda: cp1_completeness(sub, summary, report)),
        ("cp2_units", lambda: cp2_units(sub, summary, report)),
        ("cp3_common_vector", lambda: cp3_common(sub, summary)),
        ("cp4_selective_vectors", lambda: cp4_selective(sub, summary)),
        ("cp5_boundary", lambda: cp5_boundary(sub, summary, report)),
        ("cp6_bias_anomaly", lambda: cp6_bias_anomaly(sub, summary, report)),
        ("cp7_width_fourier", lambda: cp7_width_fourier(sub, summary)),
        ("cp8_spectroscopy", lambda: cp8_spectroscopy(sub, summary)),
        ("cp9_retrace", lambda: cp9_retrace(sub, summary)),
        ("cp10_report", lambda: cp10_report(sub, report)),
    ]:
        try:
            earned, mx, detail = fn()
        except Exception as exc:  # never crash the verifier
            earned, mx, detail = 0.0, 0.0, {"error": repr(exc)}
        checkpoints[key] = {"earned": round(float(earned), 3), "max": round(float(mx), 3), "detail": detail}

    total = sum(c["earned"] for c in checkpoints.values())
    max_total = sum(c["max"] for c in checkpoints.values()) or 100.0
    total = float(min(total, max_total))
    normalized = round(total / max_total, 4)

    reward_json = {"reward": normalized, "score_out_of_100": round(total, 3), "max_score": round(max_total, 3)}
    for k, c in checkpoints.items():
        reward_json[k] = c["earned"]

    (logs / "reward.json").write_text(json.dumps(reward_json, indent=2), encoding="utf-8")
    (logs / "reward.txt").write_text(f"{normalized}\n", encoding="utf-8")
    (logs / "grading_report.json").write_text(
        json.dumps({"total": round(total, 3), "max": round(max_total, 3),
                    "normalized": normalized, "checkpoints": checkpoints}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(f"[grader] submission: {sub}")
    for k, c in checkpoints.items():
        print(f"  {k:24s} {c['earned']:6.2f} / {c['max']:5.1f}")
    print(f"[grader] TOTAL {total:.2f} / {max_total:.1f}  ->  reward={normalized}")
    return normalized


def main(argv=None):
    ap = argparse.ArgumentParser(description="Harbor verifier for the GdTe2 boundary-charge-order task")
    ap.add_argument("--submission-dir", default=os.environ.get("SUBMISSION_DIR", "/app/output"))
    ap.add_argument("--logs-dir", default=os.environ.get("LOGS_DIR", "/logs/verifier"))
    args = ap.parse_args(argv)
    grade(Path(args.submission_dir), Path(args.logs_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
