"""Objective validation of the two-temperature reference results."""

from __future__ import annotations

import argparse
import ast
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _runtime_import_audit() -> tuple[bool, list[str]]:
    violations: list[str] = []
    for path in sorted((ROOT / "src").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            else:
                continue
            for name in names:
                if name == "toymodel" or name.startswith("toymodel."):
                    violations.append(f"{path.relative_to(ROOT)}:{node.lineno}:{name}")
    return not violations, violations


def _finite_weights(run: dict) -> bool:
    return all(
        math.isfinite(value)
        for chunk in run["rpsh_ba"]["chunks"]
        for value in chunk["weights"]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reference-143",
        type=Path,
        default=ROOT / "output" / "reference_143K.json",
    )
    parser.add_argument(
        "--reference-250",
        type=Path,
        default=ROOT / "output" / "reference_250K.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "output" / "validation.json",
    )
    args = parser.parse_args()

    run143 = _load(args.reference_143)
    run250 = _load(args.reference_250)
    targets = {
        143: {"exact_a": -1.1655442, "exact_na": -1.2976844, "ba": -1.30},
        250: {"exact_a": 0.1926060, "exact_na": 0.1128989, "ba": 0.14},
    }
    runs = {143: run143, 250: run250}
    import_ok, violations = _runtime_import_audit()

    checks: dict[str, bool] = {"no_external_project_import": import_ok}
    metrics: dict[str, float | int | str] = {}
    for temperature, run in runs.items():
        exact_a = run["exact"]["adiabatic"]["log10_rate_m_per_molecule_s"]
        exact_na = run["exact"]["nonadiabatic"]["log10_rate_m_per_molecule_s"]
        ba = run["rpsh_ba"]["total"]["log10_rate_m_per_molecule_s"]
        ba_se = run["rpsh_ba"]["total"]["log10_standard_error_dex"]
        key = str(temperature)
        checks[f"{key}_complete"] = run["status"] == "complete"
        checks[f"{key}_temperature"] = float(run["temperature_K"]) == float(temperature)
        checks[f"{key}_beads"] = int(run["nbeads"]) == 32
        checks[f"{key}_exact_adiabatic"] = abs(exact_a - targets[temperature]["exact_a"]) <= 0.03
        checks[f"{key}_exact_nonadiabatic"] = abs(exact_na - targets[temperature]["exact_na"]) <= 0.03
        checks[f"{key}_ba_matches_figure"] = abs(ba - targets[temperature]["ba"]) <= 0.15
        checks[f"{key}_ba_se"] = ba_se <= 0.15
        checks[f"{key}_finite_weights"] = _finite_weights(run)
        metrics[f"{key}_exact_a_log10"] = exact_a
        metrics[f"{key}_exact_na_log10"] = exact_na
        metrics[f"{key}_ba_log10"] = ba
        metrics[f"{key}_ba_se_dex"] = ba_se
        metrics[f"{key}_runtime_seconds"] = run["runtime"]["total_seconds"]
        metrics[f"{key}_python"] = run["environment"]["python"].splitlines()[0]
        metrics[f"{key}_numpy"] = run["environment"]["numpy"]
        metrics[f"{key}_scipy"] = run["environment"]["scipy"]

    checks["combined_runtime_below_120_minutes"] = (
        run143["runtime"]["total_seconds"] + run250["runtime"]["total_seconds"] < 7200.0
    )
    payload = {
        "status": "pass" if all(checks.values()) else "fail",
        "checks": checks,
        "metrics": metrics,
        "import_violations": violations,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

