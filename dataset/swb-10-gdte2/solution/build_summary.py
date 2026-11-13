#!/usr/bin/env python3
"""Map the reference pipeline's reference_results.json into the agent-facing
summary.json schema, and rename the six *_reference.png figures to the required
deliverable names.

This is a pure post-processing wrapper: every number it emits comes from
reference_results.json, which reference_analysis.py recomputed from /app/data.
Nothing is hard-coded here.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def r(x, nd):
    """Round to nd decimals if numeric, else pass through."""
    if x is None:
        return None
    try:
        return round(float(x), nd)
    except (TypeError, ValueError):
        return x


def vec(d):
    """Reference vector dict -> agent {qx, qy, magnitude} (q: 3 decimals)."""
    return {
        "qx": r(d.get("qx_cycles_per_nm"), 3),
        "qy": r(d.get("qy_cycles_per_nm"), 3),
        "magnitude": r(d.get("magnitude_cycles_per_nm"), 3),
    }


def vec_list(items):
    return [vec(d) for d in (items or [])]


# reference PNG name -> required deliverable name
PNG_RENAMES = {
    "reciprocal_space_reference.png": "reciprocal_space.png",
    "domain_segmentation_reference.png": "domain_segmentation.png",
    "bias_dependence_reference.png": "bias_dependence.png",
    "boundary_feature_fft_reference.png": "boundary_feature_fft.png",
    "spectroscopy_reference.png": "spectroscopy_comparison.png",
    "validation_reference.png": "validation_retrace.png",
}


def build_summary(ref: dict) -> dict:
    modv = ref.get("modulation_vectors", {})
    common = modv.get("common", {})
    boundary = ref.get("boundary", {})
    bias = ref.get("bias_dependence", {})
    feat = ref.get("boundary_fourier_feature", {})
    spec = ref.get("spectroscopy", {})
    retrace = ref.get("independent_retrace_validation", {})

    maps = bias.get("maps", [])
    per_bias_z = {f"{m.get('bias_mV'):.1f}": r(m.get("boundary_excess_z"), 2) for m in maps}
    active_fwhm = {
        f"{m.get('bias_mV'):.1f}": r(m.get("spatial_fwhm_nm"), 2)
        for m in maps
        if m.get("detected")
    }
    active_biases = bias.get("active_sampled_biases_mV", [])

    summary = {
        "schema_version": "1.0",
        "reciprocal_convention": "cycles/nm",
        "common_modulation_vector": vec(common),
        "left_domain_vectors": vec_list(modv.get("left_domain_selective")),
        "right_domain_vectors": vec_list(modv.get("right_domain_selective")),
        "boundary": {
            "sample_y_nm": [r(v, 2) for v in boundary.get("sample_y_nm", [])],
            "sample_x_nm": [r(v, 2) for v in boundary.get("sample_x_nm", [])],
            "uncertainty_nm": r(boundary.get("fit_uncertainty_nm"), 2),
        },
        "boundary_anomaly": {
            "detected": bool(len(active_biases) > 0),
            "active_sampled_biases_mV": [r(b, 1) for b in active_biases],
            "per_bias_boundary_excess_z": per_bias_z,
            "active_spatial_fwhm_nm": active_fwhm,
        },
        "fourier_feature": {
            "q_u_cycles_per_nm": r(feat.get("feature_q_u_cycles_per_nm"), 3),
            "q_y_cycles_per_nm": r(feat.get("feature_q_y_cycles_per_nm"), 3),
            "magnitude_cycles_per_nm": r(feat.get("feature_q_magnitude_cycles_per_nm"), 3),
            "angle_deg_from_u": r(feat.get("feature_angle_deg_from_u"), 2),
        },
        "spectroscopy": {
            "anomalous_point_id": spec.get("selected_point"),
            "feature_peak_mV": r(spec.get("feature_peak_mV"), 1),
            "left_edge_mV": r(spec.get("feature_left_edge_mV"), 1),
            "right_edge_mV": r(spec.get("feature_right_edge_mV"), 1),
            "width_mV": r(spec.get("feature_width_mV"), 1),
        },
        "retrace_validation": {
            "forward_active_biases_mV": [r(b, 1) for b in retrace.get("forward_detected_biases_mV", [])],
            "retrace_active_biases_mV": [r(b, 1) for b in retrace.get("retrace_detected_biases_mV", [])],
            "boundary_forward_retrace_rmse_nm": r(retrace.get("topography_boundary_forward_retrace_rmse_nm"), 2),
        },
    }
    return summary


def rename_pngs(out_dir: Path) -> dict:
    result = {}
    for src_name, dst_name in PNG_RENAMES.items():
        src = out_dir / src_name
        dst = out_dir / dst_name
        if src.exists():
            os.replace(src, dst)
            result[dst_name] = True
        else:
            result[dst_name] = dst.exists()
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build agent summary.json from reference_results.json")
    ap.add_argument("--reference-json", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ref = json.loads(Path(args.reference_json).read_text(encoding="utf-8"))
    summary = build_summary(ref)
    (out_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    renamed = rename_pngs(out_dir)
    print(json.dumps({"summary_json": "written", "pngs": renamed}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
