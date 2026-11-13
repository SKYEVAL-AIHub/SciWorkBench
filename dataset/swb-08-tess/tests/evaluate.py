"""Objective binary scorer for AST-TESS-001."""

import argparse
import csv
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oracle import bool_series, centroid_expectations, dilution_expectation, input_hashes_valid, leave_one_out, sha256_file, weighted_ephemeris


REQUIRED_FILES = {
    "analysis_pipeline.py", "analysis_report.md", "centroid_audit.csv", "diagnostic_summary.png",
    "event_metrics.csv", "injection_recovery.csv", "period_search.csv", "replay_manifest.json",
    "results.json", "sensitivity.csv",
}
REPLAY_STABLE_CSV = (
    "event_metrics.csv",
    "period_search.csv",
    "centroid_audit.csv",
    "injection_recovery.csv",
    "sensitivity.csv",
)
_REPLAY_CACHE = {}


def _json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None


def _csv(path):
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def _finite(frame, columns):
    if frame is None or not set(columns).issubset(frame.columns):
        return False
    try:
        return bool(np.isfinite(frame[list(columns)].to_numpy(dtype=float)).all())
    except Exception:
        return False


def _png_dimensions(path):
    try:
        with Path(path).open("rb") as handle:
            header = handle.read(24)
        if header[:8] != b"\x89PNG\r\n\x1a\n":
            return 0, 0
        return struct.unpack(">II", header[16:24])
    except Exception:
        return 0, 0


def _close(a, b, tolerance):
    try:
        return bool(np.isfinite(float(a)) and np.isfinite(float(b)) and abs(float(a) - float(b)) <= tolerance)
    except Exception:
        return False


def _expected_ids(data_dir, folder, suffix):
    return {path.name[:-len(suffix)] for path in (Path(data_dir) / folder).glob(f"*{suffix}")}


def _stable_results(value):
    if not isinstance(value, dict):
        return None
    stable = dict(value)
    stable.pop("runtime_seconds", None)
    return stable


def _frames_match(left, right):
    if left is None or right is None or list(left.columns) != list(right.columns) or len(left) != len(right):
        return False
    try:
        pd.testing.assert_frame_equal(
            left.reset_index(drop=True),
            right.reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            rtol=1e-10,
            atol=1e-12,
        )
        return True
    except AssertionError:
        return False


def _isolated_replay(data_dir, pipeline_path):
    pipeline_path = Path(pipeline_path)
    data_dir = Path(data_dir)
    if not pipeline_path.is_file() or not (data_dir / "materials_manifest.csv").is_file():
        return {"ok": False, "reason": "missing pipeline or public data manifest"}
    key = (sha256_file(pipeline_path), sha256_file(data_dir / "materials_manifest.csv"))
    if key in _REPLAY_CACHE:
        return _REPLAY_CACHE[key]

    result = {"ok": False, "reason": "isolated replay did not complete"}
    try:
        with tempfile.TemporaryDirectory(prefix="ast-tess-r26-") as tempdir:
            root = Path(tempdir)
            isolated_data = root / "data"
            isolated_pipeline = root / "analysis_pipeline.py"
            isolated_output = root / "output"
            shutil.copytree(data_dir, isolated_data)
            shutil.copy2(pipeline_path, isolated_pipeline)
            subprocess_identity = {}
            if hasattr(os, "geteuid") and os.geteuid() == 0:
                root.chmod(0o755)
                isolated_pipeline.chmod(0o644)
                for path in (isolated_data, *isolated_data.rglob("*")):
                    path.chmod(0o755 if path.is_dir() else 0o644)
                for directory in (isolated_output, root / "mplconfig"):
                    directory.mkdir()
                    os.chown(directory, 65534, 65534)
                subprocess_identity = {"user": 65534, "group": 65534, "extra_groups": []}
            environment = {"PATH": os.defpath + ":/usr/local/bin", "LANG": "C.UTF-8"}
            environment["MPLCONFIGDIR"] = str(root / "mplconfig")
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            environment["PYTHONNOUSERSITE"] = "1"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(isolated_pipeline),
                    "--data",
                    str(isolated_data),
                    "--output",
                    str(isolated_output),
                ],
                cwd=root,
                **subprocess_identity,
                env=environment,
                capture_output=True,
                text=True,
                timeout=7200,
                check=False,
            )
            if completed.returncode != 0:
                result = {"ok": False, "reason": f"pipeline exit code {completed.returncode}"}
            elif not isolated_output.is_dir():
                result = {"ok": False, "reason": "pipeline created no output directory"}
            else:
                present = {path.name for path in isolated_output.iterdir() if path.is_file()}
                replay = _json(isolated_output / "replay_manifest.json") or {}
                replay_hashes = replay.get("output_sha256", {}) if isinstance(replay, dict) else {}
                hash_names = REQUIRED_FILES - {"replay_manifest.json"}
                hashes_valid = set(replay_hashes) == hash_names and all(
                    (isolated_output / name).is_file()
                    and sha256_file(isolated_output / name) == replay_hashes[name]
                    for name in hash_names
                )
                if present != REQUIRED_FILES:
                    result = {"ok": False, "reason": f"wrong output set: {sorted(present)}"}
                elif not hashes_valid:
                    result = {"ok": False, "reason": "replayed output hashes are invalid"}
                elif not input_hashes_valid(isolated_data):
                    result = {"ok": False, "reason": "pipeline modified isolated public data"}
                else:
                    result = {
                        "ok": True,
                        "reason": "isolated replay completed",
                        "results": _stable_results(_json(isolated_output / "results.json")),
                        "tables": {name: _csv(isolated_output / name) for name in REPLAY_STABLE_CSV},
                    }
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        result = {"ok": False, "reason": f"isolated replay error: {type(error).__name__}"}
    _REPLAY_CACHE[key] = result
    return result


def _replay_matches_submission(data_dir, submission, submitted_results):
    replayed = _isolated_replay(data_dir, Path(submission) / "analysis_pipeline.py")
    if not replayed.get("ok"):
        return False
    if replayed.get("results") != _stable_results(submitted_results):
        return False
    return all(
        _frames_match(replayed["tables"].get(name), _csv(Path(submission) / name))
        for name in REPLAY_STABLE_CSV
    )


def score_submission(data_dir, submission_dir, rubric_path):
    data_dir = Path(data_dir)
    submission = Path(submission_dir)
    rubric = _json(rubric_path)
    results = _json(submission / "results.json") or {}
    replay = _json(submission / "replay_manifest.json") or {}
    events = _csv(submission / "event_metrics.csv")
    periods = _csv(submission / "period_search.csv")
    centroids = _csv(submission / "centroid_audit.csv")
    injections = _csv(submission / "injection_recovery.csv")
    sensitivity = _csv(submission / "sensitivity.csv")
    report = (submission / "analysis_report.md").read_text(encoding="utf-8", errors="ignore") if (submission / "analysis_report.md").exists() else ""
    pipeline = (submission / "analysis_pipeline.py").read_text(encoding="utf-8", errors="ignore") if (submission / "analysis_pipeline.py").exists() else ""

    present = {path.name for path in submission.iterdir()} if submission.exists() else set()
    checks = {}
    checks["R01"] = present == REQUIRED_FILES and all((submission / name).is_file() for name in REQUIRED_FILES)
    required_result_keys = {"task_id", "random_seed", "time_system", "selected_period_days", "classification", "event_count", "accepted_event_count", "centroid_event_count", "injection_count", "injection_recovery_fraction", "catalog_dilution_fraction_42arcsec", "linear_ephemeris", "leave_one_out_periods", "core_chain", "runtime_seconds"}
    checks["R02"] = required_result_keys.issubset(results) and results.get("task_id") == "AST-TESS-001" and results.get("classification") in {"robust", "ambiguous", "artifact"} and set(results.get("core_chain", {})) == {"R12", "R13", "R14", "R15", "R16", "R17"}
    checks["R03"] = input_hashes_valid(data_dir)
    replay_hashes = replay.get("output_sha256", {}) if isinstance(replay, dict) else {}
    expected_hash_names = REQUIRED_FILES - {"replay_manifest.json"}
    checks["R04"] = set(replay_hashes) == expected_hash_names and all((submission / name).is_file() and sha256_file(submission / name) == digest for name, digest in replay_hashes.items()) and all(key in replay for key in ("command", "random_seed", "bootstrap_replicates", "python", "platform", "wall_clock_seconds"))
    lower_pipeline = pipeline.lower()
    forbidden = ["/users/", "golden_output", "answer/", "grading_tool", "http://", "https://"]
    checks["R05"] = len(pipeline) > 5000 and "data_dir" in pipeline and "run_reference_analysis" in pipeline and not any(token in lower_pipeline for token in forbidden)

    event_columns = {"event_id", "epoch", "sector", "midtime_btjd", "midtime_bjd_tdb", "midtime_error_days", "depth", "depth_error", "duration_days", "exposure_seconds", "oc_seconds", "accepted"}
    event_ids = _expected_ids(data_dir, "event_windows", ".csv.gz")
    checks["R06"] = events is not None and event_columns.issubset(events.columns) and len(events) == 5 and events["event_id"].is_unique and set(events["event_id"]) == event_ids
    period_columns = {"alias_label", "period_days", "bls_snr", "bls_t0_btjd", "bls_duration_days", "bls_depth", "window_function_power", "covered_events", "supported_events", "contradictions", "event_consistency", "adjudication_score", "selected"}
    period_selected = bool_series(periods["selected"]) if periods is not None and "selected" in periods else pd.Series(dtype=bool)
    checks["R07"] = periods is not None and period_columns.issubset(periods.columns) and len(periods) == 3 and set(periods["alias_label"]) == {"P/2", "P", "2P"} and period_selected.fillna(False).sum() == 1
    centroid_columns = {"event_id", "difference_x_pix", "difference_y_pix", "reference_x_pix", "reference_y_pix", "offset_pix", "offset_arcsec", "difference_snr", "catalog_dilution_fraction_42arcsec", "on_target"}
    checks["R08"] = centroids is not None and centroid_columns.issubset(centroids.columns) and len(centroids) == 5 and centroids["event_id"].is_unique and checks["R06"] and set(centroids["event_id"]) == set(events["event_id"]) and _finite(centroids, ["offset_pix", "offset_arcsec"]) and (centroids[["offset_pix", "offset_arcsec"]].to_numpy() >= 0).all()
    injection_columns = {"injection_id", "depth", "duration_days", "period_days", "phase_offset_days", "recovered_period_days", "recovered_t0_btjd", "recovered_depth", "detection_snr", "phase_error_days", "recovered"}
    public_grid = _csv(data_dir / "config" / "injection_grid.csv")
    checks["R09"] = injections is not None and injection_columns.issubset(injections.columns) and len(injections) == 120 and injections["injection_id"].is_unique and public_grid is not None and set(injections["injection_id"]) == set(public_grid["injection_id"])
    if checks["R09"]:
        joined = injections.merge(public_grid, on="injection_id", suffixes=("", "_expected"))
        checks["R09"] = all(np.allclose(joined[name], joined[f"{name}_expected"], rtol=0, atol=1e-12) for name in ("depth", "duration_days", "period_days", "phase_offset_days"))
    scenario_columns = {"scenario_id", "quality_bits", "flare_sigma", "detrend_days", "exposure_integrated", "selected_period_days", "selected_alias", "supported_events", "classification_stable"}
    public_scenarios = _csv(data_dir / "config" / "sensitivity_scenarios.csv")
    checks["R10"] = sensitivity is not None and scenario_columns.issubset(sensitivity.columns) and len(sensitivity) == 8 and sensitivity["scenario_id"].is_unique and public_scenarios is not None and set(sensitivity["scenario_id"]) == set(public_scenarios["scenario_id"])
    if checks["R10"]:
        joined = sensitivity.merge(public_scenarios, on="scenario_id", suffixes=("", "_expected"))
        checks["R10"] = all(np.allclose(joined[name], joined[f"{name}_expected"], rtol=0, atol=1e-12) for name in ("quality_bits", "flare_sigma", "detrend_days"))
    width, height = _png_dimensions(submission / "diagnostic_summary.png")
    checks["R11"] = width >= 1200 and height >= 700 and (submission / "diagnostic_summary.png").stat().st_size > 40000 if (submission / "diagnostic_summary.png").exists() else False

    r12 = results.get("core_chain", {}).get("R12", {})
    checks["R12"] = checks["R06"] and "BTJD" in str(results.get("time_system")) and "BJD_TDB" in str(results.get("time_system")) and r12.get("quality_reject_bits") == 175 and 30000 <= int(r12.get("accepted_rows", 0)) <= 40000 and 20 <= int(r12.get("flare_rows", 0)) <= 3000 and _finite(events, ["midtime_btjd", "midtime_bjd_tdb", "exposure_seconds"]) and np.allclose(events["midtime_bjd_tdb"] - events["midtime_btjd"], 2457000.0, rtol=0, atol=1e-8) and np.allclose(events["exposure_seconds"], 120.0, rtol=0, atol=1e-6)
    r13 = results.get("core_chain", {}).get("R13", {})
    accepted_series = bool_series(events["accepted"]) if events is not None and "accepted" in events else pd.Series(dtype=bool)
    raw_r13 = checks["R06"] and all(str(r13.get(key, "")).strip() for key in ("method_primary", "method_secondary")) and _finite(events, ["depth", "depth_error", "duration_days"]) and bool(((events["depth"] > 0) & (events["depth"] < 0.02) & (events["depth_error"] > 0) & events["duration_days"].between(0.05, 0.20)).all()) and accepted_series.fillna(False).sum() >= 4
    checks["R13"] = checks["R12"] and raw_r13
    seed = _json(data_dir / "config" / "event_seed_schedule.json") or {}
    raw_r14 = checks["R07"] and _finite(periods, ["period_days", "bls_snr", "window_function_power", "supported_events", "contradictions", "event_consistency", "adjudication_score"]) and periods["window_function_power"].between(0, 1).all()
    if raw_r14:
        expected_scores = 3 * periods["supported_events"] - 2 * periods["contradictions"] + 2 * periods["event_consistency"]
        selected_row = periods.loc[period_selected.fillna(False)].iloc[0]
        raw_r14 = np.allclose(expected_scores, periods["adjudication_score"], atol=1e-9, rtol=0) and selected_row["alias_label"] == "P" and abs(float(selected_row["period_days"]) - float(seed.get("period_seed_days", -99))) < 0.03 and int(np.argmax(periods["adjudication_score"].to_numpy())) == int(np.where(period_selected.fillna(False).to_numpy())[0][0])
    checks["R14"] = checks["R13"] and bool(raw_r14)
    raw_r15 = checks["R06"] and _finite(events, ["epoch", "midtime_btjd", "midtime_error_days", "oc_seconds"]) and (events["midtime_error_days"] > 0).all() and int(results.get("core_chain", {}).get("R15", {}).get("bootstrap_replicates", 0)) >= 300
    if raw_r15:
        expected_fit = weighted_ephemeris(events)
        submitted_fit = results.get("linear_ephemeris", {})
        expected_loo = leave_one_out(events)
        submitted_loo = results.get("leave_one_out_periods", [])
        public_events = {item["event_id"]: item["center_btjd"] for item in seed.get("events", [])}
        within_windows = all(event_id in public_events and abs(float(midtime) - float(public_events[event_id])) <= 0.05 for event_id, midtime in zip(events["event_id"], events["midtime_btjd"]))
        raw_r15 = within_windows and _close(expected_fit["t0_btjd"], submitted_fit.get("t0_btjd"), 5e-7) and _close(expected_fit["period_days"], submitted_fit.get("period_days"), 5e-7) and np.allclose(expected_fit["oc_seconds"], events["oc_seconds"], atol=0.2, rtol=0) and len(submitted_loo) == 5 and all(a["excluded_event_id"] == b["excluded_event_id"] and _close(a["period_days"], b.get("period_days"), 5e-7) for a, b in zip(expected_loo, submitted_loo))
    checks["R15"] = checks["R14"] and bool(raw_r15)
    raw_r16 = checks["R08"]
    if raw_r16:
        expected_centroids = centroid_expectations(data_dir)
        for row in centroids.itertuples(index=False):
            expected = expected_centroids.get(row.event_id, {})
            raw_r16 = raw_r16 and all(_close(getattr(row, key), expected.get(key), 0.10 if "pix" in key else 2.1) for key in ("difference_x_pix", "difference_y_pix", "reference_x_pix", "reference_y_pix", "offset_pix", "offset_arcsec"))
        expected_dilution = dilution_expectation(data_dir)
        raw_r16 = raw_r16 and np.allclose(centroids["catalog_dilution_fraction_42arcsec"], expected_dilution, atol=1e-6, rtol=0) and bool_series(centroids["on_target"]).fillna(False).sum() >= 4 and _close(results.get("catalog_dilution_fraction_42arcsec"), expected_dilution, 1e-6)
    checks["R16"] = checks["R15"] and bool(raw_r16)
    raw_r17 = checks["R09"] and checks["R10"] and _finite(injections, ["recovered_period_days", "phase_error_days", "duration_days", "detection_snr"])
    if raw_r17:
        expected_recovered = (np.abs(injections["recovered_period_days"] - float(seed["period_seed_days"])) < 0.03) & (injections["phase_error_days"] <= injections["duration_days"] / 2) & (injections["detection_snr"] >= 7.0)
        actual_recovered = bool_series(injections["recovered"]).fillna(False)
        stable = bool_series(sensitivity["classification_stable"]).fillna(False)
        recovery_fraction = float(actual_recovered.mean())
        raw_r17 = bool((expected_recovered == actual_recovered).all()) and stable.sum() >= 6 and _close(results.get("injection_recovery_fraction"), recovery_fraction, 1e-12) and results.get("classification") == ("robust" if recovery_fraction >= 0.55 and stable.sum() >= 6 else "ambiguous")
    checks["R17"] = checks["R16"] and bool(raw_r17)

    lower_report = report.lower()
    report_terms = (("限制", "limitation"), ("盒", "box"), ("相关噪声", "correlated noise"))
    checks["R18"] = results.get("classification", "__") in lower_report and all(
        term in report for term in ("R12", "R13", "R14", "R15", "R16", "R17")
    ) and all(any(term in lower_report for term in alternatives) for alternatives in report_terms)
    checks["R19"] = checks["R06"] and checks["R08"] and checks["R09"] and results.get("event_count") == len(events) and results.get("accepted_event_count") == int(accepted_series.fillna(False).sum()) and results.get("centroid_event_count") == len(centroids) and results.get("injection_count") == len(injections)
    checks["R20"] = checks["R07"] and _close(results.get("selected_period_days"), periods.loc[period_selected.fillna(False), "period_days"].iloc[0], 1e-12) and results.get("core_chain", {}).get("R14", {}).get("selected_alias") == periods.loc[period_selected.fillna(False), "alias_label"].iloc[0]
    checks["R21"] = bool(raw_r15) and "day" in report.lower() and ("second" in report.lower() or "秒" in report)
    r16 = results.get("core_chain", {}).get("R16", {})
    checks["R22"] = checks["R08"] and _close(r16.get("median_offset_pix"), centroids["offset_pix"].median(), 1e-12) and r16.get("on_target_events") == int(bool_series(centroids["on_target"]).fillna(False).sum()) and np.allclose(centroids["catalog_dilution_fraction_42arcsec"], results.get("catalog_dilution_fraction_42arcsec"), atol=1e-12, rtol=0)
    checks["R23"] = checks["R09"] and len(set(zip(injections["depth"], injections["duration_days"], injections["phase_offset_days"]))) == 120 and _close(results.get("core_chain", {}).get("R17", {}).get("recovery_fraction"), bool_series(injections["recovered"]).fillna(False).mean(), 1e-12)
    checks["R24"] = checks["R10"] and results.get("core_chain", {}).get("R17", {}).get("stable_scenarios") == int(bool_series(sensitivity["classification_stable"]).fillna(False).sum()) and {0, 175}.issubset(set(sensitivity["quality_bits"])) and {False, True}.issubset(set(bool_series(sensitivity["exposure_integrated"]).dropna()))
    checks["R25"] = all(term in lower_report for term in ("mast", "tess", "gilbert")) and ("2109.03924" in report or "10.3847/1538-3881/ac23ca" in lower_report)
    checks["R26"] = (
        0 < float(results.get("runtime_seconds", 0)) <= 7200
        and 0 < float(replay.get("wall_clock_seconds", 0)) <= 7200
        and checks["R04"]
        and checks["R05"]
        and _replay_matches_submission(data_dir, submission, results)
    )

    triggered_gates = []
    if not checks["R03"]:
        triggered_gates.append("G01")
    unsafe_tokens = ["pickle.loads", "yaml.load(", "os.system(", "eval(", "exec("]
    if any(token in lower_pipeline for token in unsafe_tokens):
        triggered_gates.append("G02")
    if any(token in lower_pipeline for token in ("/users/", "golden_output", "answer/", "grading_tool")):
        triggered_gates.append("G03")
    if len(pipeline) < 1000:
        triggered_gates.append("G04")

    score_items = []
    score_raw = 0
    for item in rubric["items"]:
        if item.get("gating"):
            passed = item["id"] not in triggered_gates
            score_items.append({"id": item["id"], "passed": passed, "awarded": 0})
        else:
            passed = bool(checks.get(item["id"], False))
            awarded = item["score"] if passed else 0
            score_raw += awarded
            score_items.append({"id": item["id"], "passed": passed, "awarded": awarded})
    return {
        "task_id": "AST-TESS-001",
        "score_raw": int(score_raw),
        "score_norm": score_raw / 100.0,
        "score_final": 0.0 if triggered_gates else score_raw / 100.0,
        "triggered_gates": triggered_gates,
        "items": score_items,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--submission", type=Path, required=True)
    parser.add_argument("--rubric", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    result = score_submission(args.data, args.submission, args.rubric)
    text = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.json_out:
        args.json_out.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
