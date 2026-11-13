"""Compare guide-level tables to the frozen reference estimates."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd


def read(path, keys):
    data = pd.read_csv(path, sep="\t")
    data = data.set_index(list(keys))
    if not data.index.is_unique:
        raise ValueError("Duplicate row keys")
    return data


def similarity(submission, expected, keys, columns):
    try:
        actual = read(submission, keys).reindex(expected.index)
    except (OSError, ValueError, KeyError):
        return 0.0
    score = 0.0
    for column, tolerance in columns.items():
        if column not in actual:
            continue
        a, b = actual[column], expected[column]
        present = b.notna()
        if not present.any():
            continue
        a, b = a[present], b[present]
        if tolerance is None:
            score += (a.fillna("").astype(str).str.lower() == b.astype(str).str.lower()).mean()
        else:
            a = pd.to_numeric(a, errors="coerce")
            if column in {"PC1", "PC2"}:
                a = -a if (a * b).sum() < 0 else a
            score += (np.isfinite(a) & ((a - b).abs() <= tolerance)).mean()
    return score / len(columns)


def grade(output, reference):
    checks = [
        ("ago2_gene_effects.tsv.gz", ("gene_symbol",), {
            "eligible_main": None, "eligible_strict": None, "AGO2_guide_id": None,
            "whole_log2_effect": .15, "nascent_log2_effect": .2, "preexisting_log2_effect": .15,
            "common_scale_total_effect": .15, "nascent_shapley_contribution": .15,
            "preexisting_shapley_contribution": .15, "decomposition_residual": .05,
        }, 40),
        ("perturbation_comparison.tsv", ("target_gene",), {
            "n_guides": 0, "rho_whole_nascent": .1, "rho_whole_preexisting": .1,
            "nascent_fraction_ge_0_5": .07, "median_target_whole_log2_effect": .3,
        }, 14),
        ("placebo_calibration.tsv", ("analysis_level", "query_id"), {
            "reference_profiles": 0, "rho_whole_nascent": .1, "rho_whole_preexisting": .1,
            "nascent_fraction_ge_0_5": .07, "query_class": None,
        }, 21),
        ("guide_qc.tsv", ("guide_id",), {
            "target_gene": None, "cells_aggregated": 0, "low_cell_count_flag": None,
            "nascent_fraction_of_whole": .02, "whole_library_size": 1,
        }, 8),
        ("pca_coordinates.tsv", ("guide_id",), {
            "PC1_variance_fraction": .08, "PC2_variance_fraction": .06, "PC1": 2, "PC2": 2,
        }, 6),
        ("sensitivity_analysis.tsv", ("configuration",), {
            "eligible_genes": 150, "whole_effect_spearman_vs_main": .12,
            "nascent_effect_spearman_vs_main": .12, "preexisting_effect_spearman_vs_main": .12,
        }, 7),
    ]
    score = 0.0
    for name, keys, columns, weight in checks:
        expected = read(reference / f"expected_{name}", keys)
        if name.startswith("ago2_gene_effects"):
            expected = expected[expected["eligible_main"].astype(bool)]
        score += weight * similarity(output / name, expected, keys, columns)
    score += 4 if (output / "claim_evidence.tsv").is_file() and (output / "report.md").is_file() else 0
    return float(score / 100)


if __name__ == "__main__":
    print(grade(*(Path(arg) for arg in sys.argv[1:])))
