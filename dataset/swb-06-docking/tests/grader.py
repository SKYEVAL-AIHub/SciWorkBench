#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("/app/output"))
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--reward", type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    try:
        data = json.loads((out / "review_findings.json").read_text(encoding="utf-8"))
    except Exception:
        data = {}
    try:
        report = (out / "review_report.md").read_text(encoding="utf-8").lower()
    except Exception:
        report = ""
    checks = []

    def check(module, name, maximum, condition, actual=None):
        checks.append({"module": module, "check": name, "points": float(maximum if condition else 0), "maximum": float(maximum), "passed": bool(condition), "actual": actual})

    def near(module, name, maximum, key, target, tolerance):
        actual = data.get(key)
        try:
            ok = abs(float(actual) - target) <= tolerance
        except (TypeError, ValueError):
            ok = False
        check(module, name, maximum, ok, {"value": actual, "expected": target})

    def rmsd(module, name, maximum, key):
        value = data.get(key)
        try:
            value = float(value)
            points = maximum if value < 3.0 else (0.5 if value <= 5.0 else 0.0)
        except (TypeError, ValueError):
            points = 0.0
        checks.append({"module": module, "check": name, "points": points, "maximum": maximum, "passed": points == maximum, "actual": value})

    m = "1_independent_docking"
    rmsd(m, "CDK2 redocking RMSD", 1.5, "cdk2_redocking_rmsd")
    rmsd(m, "ROCK1 redocking RMSD", 1.5, "rock1_redocking_rmsd")
    check(m, "docking result table", 2, (out / "docking_results.csv").is_file())
    check(m, "computed properties", 1, (out / "property_results.csv").is_file())
    check(m, "critical comparison in report", 2,
          ("独立" in report or "independent" in report) and
          any(word in report for word in ("差异", "比对", "difference", "comparison", "agreement")))

    m = "2_performance_and_statistics"
    directions = data.get("metric_directions", {})
    expected_directions = {"Vina Score": "lower", "Vina Min": "lower", "Vina Dock": "lower", "QED": "higher", "SA": "higher", "High_Affinity": "higher", "Lipinski": "higher", "Diversity": "higher"}
    check(m, "metric directions", 4, directions == expected_directions, directions)
    check(m, "best method is not uniformly best", 3, data.get("method_i_best_vina") is True and data.get("method_i_not_best_sa_or_lipinski") is True)
    check(m, "significant methods", 3, set(data.get("significant_methods", [])) == {"Method_G", "Method_H", "Method_I"})
    check(m, "one-sided interpretation", 2, "单侧" in report or "one-sided" in report)
    check(m, "negative HLE intervals", 3, set(data.get("hle_negative_methods", [])) == {"Method_G", "Method_H", "Method_I"})
    check(m, "Method_H sign reversal", 5, data.get("sign_reversal_method") == "Method_H")
    near(m, "Method_H aggregate delta", 2, "method_h_aggregate_delta", -0.30, 0.005)
    near(m, "Method_H MPD", 1, "method_h_mpd", 0.042, 0.0005)
    check(m, "cross-level interpretation", 3, "method_h" in report and ("辛普森" in report or "异质" in report or "simpson" in report))

    m = "3_ablation_and_synergy"
    degradation = data.get("degradation_vina", {})
    expected_degradation = {"w/o_S": 1.42, "w/o_PQ": 1.38, "w/o_V": 1.07, "w/o_P": 0.68, "w/o_Q": 0.60}
    try:
        degradation_ok = set(degradation) == set(expected_degradation) and all(abs(float(degradation[k]) - v) <= 0.005 for k, v in expected_degradation.items())
    except (TypeError, ValueError):
        degradation_ok = False
    check(m, "component degradation and ranking", 5, degradation_ok, degradation)
    near(m, "Vina P-Q synergy", 4, "pq_synergy_vina_score", 0.10, 0.005)
    check(m, "FDR result", 5, data.get("fdr_significant_metrics") == ["Vina Score"], data.get("fdr_significant_metrics"))
    check(m, "multiple-testing correction implemented and discussed", 4, "fdr" in report and ("校正" in report or "correction" in report))
    near(m, "SA P gain", 1.5, "sa_p_gain", -0.0393, 0.00005)
    near(m, "SA Q gain", 1.5, "sa_q_gain", -0.0323, 0.00005)
    near(m, "SA joint gain and pattern", 3, "sa_both_gain", 0.0121, 0.00005)
    check(m, "optimal timestep", 2, data.get("optimal_timestep") == 50)
    check(m, "backbone independence boundary", 2, data.get("backbone_independence_complete") is False)

    m = "4_sr_validation_and_limits"
    check(m, "nine SR criteria", 3, data.get("sr_criteria_count") == 9)
    check(m, "C1 paired calibration impossible", 5, data.get("c1_pairwise_calibration_feasible") is False)
    normalized_report = report.replace(" ", "_")
    check(m, "missing-score mechanism explained", 4,
          "vina_score" in normalized_report and "vina_dock" in normalized_report and
          any(word in report for word in ("没有", "互斥", "缺失", "missing", "mutually exclusive")))
    check(m, "SR table and Fsp3 failure", 3, (out / "sr_screening.csv").is_file() and data.get("active_fsp3_systematic_failure") is True)
    check(m, "selected case-study molecules", 2, data.get("case_study_selected") is True)
    near(m, "full SR rate", 2, "full_sr_rate", 0.1153, 0.00005)
    check(m, "selection-bias interpretation", 2, "抽样偏差" in report or "选择偏差" in report or "selection bias" in report)
    near(m, "SA target", 1, "sa_target", 1.0, 1e-9)
    near(m, "SA actual", 1, "sa_actual", 0.60, 1e-9)
    near(m, "SA threshold", 1, "sa_sr_threshold", 0.59, 1e-9)

    m = "5_critical_review"
    near(m, "SA boundary risk", 3, "sa_actual", 0.60, 1e-9)
    check(m, "cross-level risk", 2, data.get("sign_reversal_method") == "Method_H")
    check(m, "FDR risk", 2, data.get("fdr_significant_metrics") == ["Vina Score"])
    near(m, "SA abnormal synergy", 1, "sa_both_gain", 0.0121, 0.00005)
    check(m, "timestep paradox", 0.5, data.get("optimal_timestep") == 50)
    check(m, "case selection risk", 0.5, data.get("case_study_selected") is True)
    check(m, "most vulnerable claim", 3, data.get("most_vulnerable_claim") in {"full_sr", "multi_metric_synergy"})
    check(m, "specific supplementary experiment", 2, len(report) > 1000 and any(
        phrase in report for phrase in ("补充实验", "supplementary experiment", "follow-up experiment", "additional experiment")))

    score = round(sum(c["points"] for c in checks), 4)
    result = {"benchmark": "molecular_generation_review_harbor_v1", "score": score, "maximum": 100, "checks": checks}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.reward.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    args.reward.write_text(f"{score / 100:.8f}\n", encoding="utf-8")
    print(json.dumps({"score": score, "maximum": 100}))


if __name__ == "__main__":
    main()
