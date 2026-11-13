#!/usr/bin/env python3
"""Deterministic grader for the Tully 2012 RPSH Fig. 2 rate-reproduction task.

Encodes the 8-checkpoint, 100-point rubric in grading_no_tool.md. Scores the
final products in the submission ``output/`` directory only (result-based
grading; tool-call history is not inspected).

Reference centres are the paper Fig. 2 target points (log10 k in
m molecule^-1 s^-1):

    143 K  exact_adiabatic    -1.166     exact_nonadiabatic  -1.298
    143 K  rpsh_ba (P=32)     -1.30
    250 K  exact_adiabatic     0.193     exact_nonadiabatic   0.113
    250 K  rpsh_ba (P=32)      0.14

Usage
-----
    python3 grader.py --submission-dir <output> --output-json <reward.json>
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import struct
import sys
from pathlib import Path

# --------------------------------------------------------------------------- #
# Rubric constants
# --------------------------------------------------------------------------- #
TARGETS = {
    (143.0, "exact_adiabatic"): -1.166,
    (143.0, "exact_nonadiabatic"): -1.298,
    (250.0, "exact_adiabatic"): 0.193,
    (250.0, "exact_nonadiabatic"): 0.113,
}
BA_TARGETS = {(143.0, "rpsh_ba"): -1.30, (250.0, "rpsh_ba"): 0.14}

REQUIRED_COLUMNS = [
    "temperature_K", "method", "beads", "rate_au",
    "rate_m_per_molecule_s", "log10_rate_m_per_molecule_s",
    "se_log10_dex", "numerically_converged",
]
MANDATORY_FILES = [
    "results.csv", "figure2_subset.png", "method.md",
    "conclusions.md", "diagnostics.json", "run_all.py",
]

DIAG_TOP = ["model_checks", "exact_quantum_checks", "rpsh_ba_checks", "global_checks"]
DIAG_EXACT = ["per_temperature", "spatial_boundary_check", "energy_quadrature_check",
              "open_channel_flux_check", "thermal_tail_check"]
DIAG_RPSH = ["per_temperature", "sampling_size", "independent_blocks",
             "standard_error_estimator", "reactive_trajectory_count",
             "crossing_count_summary", "hop_count_summary",
             "electronic_norm_error", "energy_error"]

PNG_SIG = b"\x89PNG\r\n\x1a\n"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _to_float(x):
    try:
        v = float(str(x).strip())
        return v
    except Exception:
        return None


def _png_dpi(path: Path):
    """Read the pHYs chunk to recover the PNG resolution in dpi (None if absent)."""
    try:
        data = path.read_bytes()
    except Exception:
        return None
    if not data.startswith(PNG_SIG):
        return None
    i = 8
    n = len(data)
    while i + 8 <= n:
        length = struct.unpack(">I", data[i:i + 4])[0]
        ctype = data[i + 4:i + 8]
        body = data[i + 8:i + 8 + length]
        if ctype == b"pHYs" and len(body) >= 9:
            ppu_x, ppu_y, unit = struct.unpack(">IIB", body[:9])
            if unit == 1:  # metres
                return round(ppu_x * 0.0254, 1)
            return None
        if ctype == b"IEND":
            break
        i += 12 + length
    return None


def _has_png(path: Path) -> bool:
    try:
        return path.read_bytes()[:8] == PNG_SIG
    except Exception:
        return False


def _finite_json(obj) -> bool:
    if isinstance(obj, float):
        return math.isfinite(obj)
    if isinstance(obj, dict):
        return all(_finite_json(v) for v in obj.values())
    if isinstance(obj, list):
        return all(_finite_json(v) for v in obj)
    return True


def _read_results(sub: Path):
    """Return (rows, header_ok, parse_error). rows keyed by (temp, method)."""
    path = sub / "results.csv"
    if not path.exists():
        return {}, False, "results.csv missing"
    try:
        with path.open(newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            fieldnames = reader.fieldnames or []
            raw = list(reader)
    except Exception as exc:  # noqa: BLE001
        return {}, False, f"results.csv parse error: {exc}"
    header_ok = all(c in fieldnames for c in REQUIRED_COLUMNS)
    rows = {}
    for r in raw:
        t = _to_float(r.get("temperature_K"))
        m = (r.get("method") or "").strip()
        if t is None:
            continue
        rows[(t, m)] = r
    return rows, header_ok, None


def _tier_exact(dev: float) -> int:
    a = abs(dev)
    if a <= 0.08:
        return 8
    if a <= 0.15:
        return 5
    if a <= 0.25:
        return 2
    return 0


def _tier_ba(dev: float) -> int:
    a = abs(dev)
    if a <= 0.15:
        return 14
    if a <= 0.25:
        return 9
    if a <= 0.40:
        return 4
    return 0


# --------------------------------------------------------------------------- #
# Checkpoints
# --------------------------------------------------------------------------- #
def cp1_structure(sub: Path, rows, header_ok, parse_err):
    """10 pts - mandatory artefacts + table structure."""
    detail = {}
    present = {f: (sub / f).exists() for f in MANDATORY_FILES}
    src_ok = (sub / "src").is_dir() and any((sub / "src").rglob("*.py"))
    png_ok = _has_png(sub / "figure2_subset.png")

    core_missing = [f for f in ("results.csv", "figure2_subset.png", "conclusions.md") if not present[f]]
    if core_missing or not src_ok or not present["run_all.py"]:
        return 0, {"missing_core": core_missing, "src_ok": src_ok,
                   "run_all_ok": present["run_all.py"]}

    # exactly 6 rows covering 2 temps x 3 methods
    expected_keys = {(t, m) for t in (143.0, 250.0)
                     for m in ("exact_adiabatic", "exact_nonadiabatic", "rpsh_ba")}
    row_ok = len(rows) == 6 and set(rows.keys()) == expected_keys
    beads_ok = True
    for (t, m), r in rows.items():
        b = (r.get("beads") or "").strip()
        if m == "rpsh_ba" and _to_float(b) != 32:
            beads_ok = False
        if m in ("exact_adiabatic", "exact_nonadiabatic") and b not in ("", "nan", "none", "None"):
            beads_ok = False

    detail = {"files_present": all(present.values()), "src_ok": src_ok,
              "png_ok": png_ok, "header_ok": header_ok, "row_ok": row_ok,
              "beads_ok": beads_ok, "parse_err": parse_err}

    if png_ok and header_ok and row_ok and beads_ok and all(present.values()):
        return 10, detail
    # one repairable format problem / missing non-core note
    if png_ok and row_ok:
        return 6, detail
    return 6 if png_ok else 0, detail


def cp_exact(rows, temp, cp_label):
    """16 pts - two exact quantum points at one temperature."""
    total = 0
    detail = {}
    for method in ("exact_adiabatic", "exact_nonadiabatic"):
        target = TARGETS[(temp, method)]
        r = rows.get((temp, method))
        val = _to_float(r.get("log10_rate_m_per_molecule_s")) if r else None
        if val is None or not math.isfinite(val):
            detail[method] = {"value": val, "target": target, "score": 0}
            continue
        sc = _tier_exact(val - target)
        total += sc
        detail[method] = {"value": round(val, 4), "target": target,
                          "dev_dex": round(val - target, 4), "score": sc}
    return total, detail


def cp_ba(rows, temp):
    """14 pts - one RPSH-BA point at a temperature."""
    target = BA_TARGETS[(temp, "rpsh_ba")]
    r = rows.get((temp, "rpsh_ba"))
    if not r:
        return 0, {"target": target, "reason": "row missing"}
    val = _to_float(r.get("log10_rate_m_per_molecule_s"))
    beads = _to_float((r.get("beads") or "").strip())
    se = _to_float(r.get("se_log10_dex"))
    rate = _to_float(r.get("rate_m_per_molecule_s"))
    conv = (r.get("numerically_converged") or "").strip().lower()

    # zero rate flagged converged -> 0
    if rate is not None and rate == 0.0 and conv in ("true", "1", "yes"):
        return 0, {"target": target, "value": val, "reason": "zero rate marked converged"}
    if beads != 32:
        return 0, {"target": target, "value": val, "beads": beads, "reason": "beads != 32"}
    if val is None or not math.isfinite(val):
        return 0, {"target": target, "value": val, "reason": "non-finite log10"}

    base = _tier_ba(val - target)
    # full 14-tier additionally requires a finite reported standard error
    if base == 14 and not (se is not None and math.isfinite(se) and se > 0.0):
        base = 9
    return base, {"target": target, "value": round(val, 4),
                  "dev_dex": round(val - target, 4), "beads": beads,
                  "se_log10_dex": se, "score": base}


def cp6_figure(sub: Path, rows):
    """10 pts - Fig. 2 subset plot."""
    fig = sub / "figure2_subset.png"
    if not _has_png(fig):
        return 0, {"reason": "figure2_subset.png missing or not a valid PNG"}
    dpi = _png_dpi(fig)
    # six finite plottable points
    finite_pts = 0
    for key in [(t, m) for t in (143.0, 250.0)
                for m in ("exact_adiabatic", "exact_nonadiabatic", "rpsh_ba")]:
        r = rows.get(key)
        if r:
            v = _to_float(r.get("log10_rate_m_per_molecule_s"))
            if v is not None and math.isfinite(v):
                finite_pts += 1
    dpi_ok = dpi is not None and dpi >= 150
    detail = {"dpi": dpi, "dpi_ok": dpi_ok, "finite_points": finite_pts}
    if dpi_ok and finite_pts == 6:
        return 10, detail
    if finite_pts == 6:
        return 6, detail  # data consistent but dpi unverifiable / < 150
    return 3, detail


def cp7_conclusions(sub: Path, rows):
    """10 pts - literature conclusions reproduced with numbers."""
    path = sub / "conclusions.md"
    if not path.exists():
        return 0, {"reason": "conclusions.md missing"}
    text = path.read_text(encoding="utf-8", errors="replace")
    low = text.lower()

    # direction: k_NA < k_A at both temperatures (nonadiabatic lowers the rate)
    dir_kw = any(k in low for k in [
        "k_{\\mathrm{na}}<k_{\\mathrm a}", "k_na<k_a", "k_na < k_a",
        "nonadiabatic", "非绝热", "降低", "小于", "lower",
    ])
    both_temps = ("143" in text) and ("250" in text)
    ba_kw = any(k in low for k in ["rpsh", "ba", "bead", "复现", "reproduce", "一致"])
    # numeric citations from results.csv
    nums = []
    for key, tgt in list(TARGETS.items()) + list(BA_TARGETS.items()):
        r = rows.get(key)
        if r:
            v = _to_float(r.get("log10_rate_m_per_molecule_s"))
            if v is not None:
                nums.append(f"{v:.3f}")
    cited = sum(1 for s in nums if s in text)
    detail = {"direction_kw": dir_kw, "both_temps": both_temps,
              "ba_kw": ba_kw, "numeric_citations": cited, "of": len(nums)}

    if dir_kw and both_temps and ba_kw and cited >= 3:
        return 10, detail
    if dir_kw and ba_kw and (both_temps or cited >= 1):
        return 6, detail
    if dir_kw or ba_kw:
        return 3, detail
    return 0, detail


def cp8_diagnostics(sub: Path):
    """10 pts - diagnostics.json + method.md numerical checks."""
    dpath = sub / "diagnostics.json"
    mpath = sub / "method.md"
    detail = {}
    if not dpath.exists():
        return 0, {"reason": "diagnostics.json missing"}
    try:
        diag = json.loads(dpath.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        return 0, {"reason": f"diagnostics.json parse error: {exc}"}

    top_ok = all(k in diag for k in DIAG_TOP)
    finite_ok = _finite_json(diag)
    ex = diag.get("exact_quantum_checks", {})
    rp = diag.get("rpsh_ba_checks", {})
    ex_ok = all(k in ex for k in DIAG_EXACT)
    rp_ok = all(k in rp for k in DIAG_RPSH)
    # RPSH per-temperature block should carry the sampling/hop/crossing detail
    per_t = rp.get("per_temperature", {})
    per_t_ok = isinstance(per_t, dict) and len(per_t) >= 2

    method_txt = mpath.read_text(encoding="utf-8", errors="replace").lower() if mpath.exists() else ""
    method_kw = sum(1 for k in ["boundary", "边界", "quadrature", "求积", "积分",
                                "sampling", "取样", "sampling size", "standard error",
                                "标准误", "open channel", "开放通道", "chunk", "分块"]
                    if k in method_txt)

    detail = {"top_ok": top_ok, "finite_ok": finite_ok, "exact_ok": ex_ok,
              "rpsh_ok": rp_ok, "per_temperature_ok": per_t_ok,
              "method_keywords": method_kw, "method_present": mpath.exists()}

    if top_ok and finite_ok and ex_ok and rp_ok and per_t_ok and method_kw >= 3:
        return 10, detail
    if top_ok and finite_ok and (ex_ok or rp_ok) and method_kw >= 1:
        return 6, detail
    if finite_ok and top_ok:
        return 3, detail
    return 0, detail


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def grade(sub: Path):
    rows, header_ok, parse_err = _read_results(sub)
    report = {}
    total = 0.0

    s, d = cp1_structure(sub, rows, header_ok, parse_err)
    report["checkpoint_1_structure"] = {"earned": s, "max": 10, "detail": d}
    total += s

    s, d = cp_exact(rows, 143.0, "cp2")
    report["checkpoint_2_exact_143K"] = {"earned": s, "max": 16, "detail": d}
    total += s

    s, d = cp_exact(rows, 250.0, "cp3")
    report["checkpoint_3_exact_250K"] = {"earned": s, "max": 16, "detail": d}
    total += s

    s, d = cp_ba(rows, 143.0)
    report["checkpoint_4_ba_143K"] = {"earned": s, "max": 14, "detail": d}
    total += s

    s, d = cp_ba(rows, 250.0)
    report["checkpoint_5_ba_250K"] = {"earned": s, "max": 14, "detail": d}
    total += s

    s, d = cp6_figure(sub, rows)
    report["checkpoint_6_figure"] = {"earned": s, "max": 10, "detail": d}
    total += s

    s, d = cp7_conclusions(sub, rows)
    report["checkpoint_7_conclusions"] = {"earned": s, "max": 10, "detail": d}
    total += s

    s, d = cp8_diagnostics(sub)
    report["checkpoint_8_diagnostics"] = {"earned": s, "max": 10, "detail": d}
    total += s

    report["total"] = round(total, 2)
    report["max_total"] = 100
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--submission-dir", required=True, type=Path)
    ap.add_argument("--output-json", type=Path, default=None)
    args = ap.parse_args()

    sub = args.submission_dir
    if not sub.is_dir():
        print(f"ERROR: submission dir not found: {sub}", file=sys.stderr)
        report = {"total": 0.0, "max_total": 100, "error": "submission dir missing"}
    else:
        report = grade(sub)

    reward = round(float(report.get("total", 0.0)) / 100.0, 4)
    payload = {"reward": reward, "score": report.get("total", 0.0),
               "max_score": 100, "report": report}

    text = json.dumps(payload, indent=2, ensure_ascii=False)
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
