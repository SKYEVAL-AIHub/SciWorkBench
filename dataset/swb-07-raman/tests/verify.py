"""Score six automated checks of the Raman submission artifacts."""

import csv
import json
import math
from pathlib import Path
import stat
import sys

from PIL import Image


def read_file(root, name):
    path = root / name
    if root.is_symlink():
        raise ValueError("Symlink output directory")
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode) or not 0 < info.st_size <= 32 * 1024 * 1024:
        raise ValueError(f"Invalid artifact: {name}")
    return path


def raman_value(output, component, metric, target, tolerance, helicity=None):
    with read_file(output, "fit_results.csv").open(newline="", encoding="utf-8-sig") as stream:
        candidates = [
            row for row in csv.DictReader(stream)
            if row["component"] == component and row["metric"] == metric
            and (helicity is None or row["helicity"] == helicity)
        ]
    if len(candidates) != 1:
        return False
    row = candidates[0]
    uncertainty = float(row["uncertainty"])
    unit = "cm^-1" if metric == "center" else "mu_B"
    return (
        row["unit"] == unit and math.isfinite(uncertainty) and uncertainty > 0
        and math.isclose(float(row["value"]), target, rel_tol=0, abs_tol=tolerance)
    )


def raman_figure(output):
    with Image.open(read_file(output, "raman_analysis.png")) as image:
        if image.format != "PNG" or min(image.size) < 300:
            return False
        if not all(295 <= value <= 305 for value in image.info.get("dpi", (0, 0))):
            return False
        image.verify()
    return True


def main():
    output, logs = map(Path, sys.argv[1:])
    original = "original_hybridized_mode"
    dopant = "dopant_induced_mode"
    checks = {
        "original_mode_center": lambda: raman_value(output, original, "center", 264.92, 1.5),
        "dopant_mode_center": lambda: raman_value(output, dopant, "center", 282.86, 1.5),
        "original_mode_moments": lambda: all(
            raman_value(output, original, "mu_eff", 2.25, 0.25, channel)
            for channel in ("RL", "LR")
        ),
        "dopant_mode_moments": lambda: all(
            raman_value(output, dopant, "mu_eff", 2.8, 0.3, channel)
            for channel in ("RL", "LR")
        ),
        "report_present": lambda: bool(read_file(output, "raman_results.md").read_text().strip()),
        "figure_present_at_requested_dpi": lambda: raman_figure(output),
    }
    results = {}
    for name, check in checks.items():
        try:
            results[name] = bool(check())
        except (OSError, ValueError, KeyError, TypeError, IndexError, csv.Error) as error:
            print(f"{name}: {type(error).__name__}")
            results[name] = False
    reward = sum(results.values()) / len(results)
    logs.mkdir(parents=True, exist_ok=True)
    (logs / "checks.json").write_text(json.dumps(results, indent=2) + "\n")
    (logs / "reward.txt").write_text(f"{reward}\n")
    print(json.dumps({"reward": reward, **results}))


if __name__ == "__main__":
    main()
