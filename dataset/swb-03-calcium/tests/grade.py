"""Score the numeric parts of the calcium-imaging rubric."""

import csv
from pathlib import Path
import sys

import numpy as np
from PIL import Image


def read_table(path, expected, keys):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(expected) or any(
        tuple(a) != tuple(b) or any(a[k] != b[k] for k in keys)
        for a, b in zip(rows, expected)
    ):
        raise ValueError("Unexpected table layout or row keys")
    return rows


def numbers(rows, columns):
    return np.array([[float(row[column] or "nan") for column in columns] for row in rows])


def within(a, b, absolute, relative=0, offset=0):
    return (np.isnan(a) & np.isnan(b)) | (
        np.isfinite(a) & np.isfinite(b) &
        (np.abs(a - b) <= np.maximum(absolute, relative * (np.abs(b) + offset)))
    )


def tier(fraction, thresholds):
    return next((points for bound, points in thresholds if fraction >= bound), 0)


def grade(output, reference):
    names = ("motion_qc.csv", "calcium_traces.csv", "summary.csv", "analysis.png")
    score = sum((output / name).is_file() for name in names)
    valid_tables = 0
    for name, keys, structure in zip(names, (("frame", "time_s"), ("frame", "time_s"), ("cell",)), (1.5, 1.5, 2)):
        with (reference / name).open(newline="") as stream:
            expected = list(csv.DictReader(stream))
        try:
            rows = read_table(output / name, expected, keys)
        except (OSError, ValueError, csv.Error):
            continue
        score += structure
        valid_tables += 1
        try:
            if name == "motion_qc.csv":
                mask = np.array([row["valid"] == "yes" for row in expected])
                cols = ("shift_x_px", "shift_y_px")
                rate = within(numbers(rows, cols)[mask], numbers(expected, cols)[mask], .6).mean()
                score += tier(rate, ((.95, 12), (.85, 9), (.70, 6), (.50, 3)))
                correct = sum(a["valid"] == b["valid"] for a, b in zip(rows, expected))
                score += tier(correct, ((238, 4), (234, 2)))
                cols = ("global_mean_ratio", "saturated_pixel_fraction")
                score += (within(numbers(rows, cols), numbers(expected, cols), .01).mean(axis=0) >= .95).sum()
            elif name == "calcium_traces.csv":
                cols = [f"cell{i}" for i in range(1, 7)]
                rate = within(numbers(rows, cols), numbers(expected, cols), 3, .05, 1).mean()
                score += tier(rate, ((.97, 32), (.92, 27), (.82, 20), (.65, 12), (.45, 6)))
            else:
                cols = ("F0", "baseline_sd_percent", "peak_dFF0_percent", "peak_time_s", "latency_half_peak_s", "AUC_percent_s")
                matches = within(numbers(rows, cols), numbers(expected, cols),
                                 [0, 1, 3, .11, .12, 0], [.04, 0, .05, 0, 0, .06], [0, 0, 1, 0, 0, 0])
                score += matches.mean(axis=0) @ [4.5, 4.5, 6, 4.5, 3, 2.5]
                score += 5 * sum(a["responder"] == b["responder"] for a, b in zip(rows, expected)) / 6
        except (TypeError, ValueError):
            continue
    score += valid_tables == 3
    try:
        with Image.open(output / "analysis.png") as image:
            if image.format == "PNG" and min(image.size) >= 300:
                image.verify()
                score += 10
    except (OSError, ValueError):
        pass
    return float(score / 100)


if __name__ == "__main__":
    print(grade(*(Path(arg) for arg in sys.argv[1:])))
