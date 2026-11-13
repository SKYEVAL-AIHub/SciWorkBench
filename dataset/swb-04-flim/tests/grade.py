"""Check coordinate coverage and validity for the sparse FLIM release."""

import csv
import math
from pathlib import Path
import sys


def coordinates(path, payload_size, header=False):
    result = {}
    with path.open(newline="") as stream:
        rows = csv.reader(line for line in stream if line.strip() and not line.lstrip().startswith("#"))
        for index, row in enumerate(rows):
            if header and index == 0 and row[0].strip().lower() in {"frame", "frame_id"}:
                continue
            if len(row) != 3 + payload_size:
                raise ValueError("Wrong column count")
            key = tuple(int(value) for value in row[:3])
            if key in result:
                raise ValueError("Duplicate coordinate")
            result[key] = tuple(
                float(value) if value.strip().lower() not in {"", "na", "n/a", "null"} else math.nan
                for value in row[3:]
            )
    return result


def valid_row(name, fields):
    # NA is allowed for derived lifetime/voltage where it cannot be estimated.
    if name == "preprocessed_flim_stack.txt":
        g1, g2, confidence = fields
        return all(math.isfinite(g) and g >= 0 for g in (g1, g2)) and confidence in (0, 1)
    if name == "tau_timeseries.txt":
        tau, error = fields
        return (math.isnan(tau) and math.isnan(error)) or (
            math.isfinite(tau) and tau > 0 and math.isfinite(error) and error >= 0
        )
    voltage, confidence, saturation = fields
    return (math.isfinite(voltage) or math.isnan(voltage)) and confidence in (0, 1) and saturation in (0, 1)


def grade(output, data):
    observed = coordinates(data / "raw_flim_stack.txt", 2)
    specs = (("preprocessed_flim_stack.txt", 3), ("tau_timeseries.txt", 2), ("Vm_mapping.txt", 3))
    score = 0.0
    for name, width in specs:
        try:
            table = coordinates(output / name, width, header=True)
        except (OSError, ValueError, IndexError, csv.Error):
            continue
        if set(table) <= set(observed):
            score += 20 * sum(valid_row(name, row) for row in table.values()) / len(observed)
    for name in ("ap_waveform.csv", "stat_ap_ratio.csv", "fig1_rmp_data.txt",
                 "fig2_ap_trace_data.txt", "fig3_apd_boxplot_data.txt", "fig4_vm_hist_data.txt"):
        path = output / name
        if path.is_file() and path.stat().st_size:
            score += 5
    report = output / "analysis_conclusion.md"
    if report.is_file():
        text = report.read_text(errors="replace").lower()
        score += 5 * any(word in text for word in ("coverage", "sparse", "missing", "insufficient"))
        score += 5 * all(word in text for word in ("calibrat", "noise", "uncertaint"))
    return score / 100


if __name__ == "__main__":
    print(grade(*(Path(arg) for arg in sys.argv[1:])))
