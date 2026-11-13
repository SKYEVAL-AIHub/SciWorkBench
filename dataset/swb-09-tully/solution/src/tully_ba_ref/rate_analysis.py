"""Combine independent transmission batches into Eq. (40) rate estimates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CombinedEstimate:
    value: float
    standard_error: float
    within_batch_standard_error: float
    between_batch_standard_error: float | None
    total_trajectories: int
    batch_count: int


def combine_independent_batches(
    values: np.ndarray, standard_errors: np.ndarray, sample_counts: np.ndarray
) -> CombinedEstimate:
    """Combine batch estimates without inverse-variance selection bias.

    Sample-count weights reconstruct the natural pooled estimator to leading
    order.  When multiple seeds are available, the larger of the propagated
    within-batch error and the seed-to-seed error is retained.
    """

    values = np.asarray(values, dtype=float)
    errors = np.asarray(standard_errors, dtype=float)
    counts = np.asarray(sample_counts, dtype=float)
    if not (values.shape == errors.shape == counts.shape) or values.ndim != 1:
        raise ValueError("values, standard_errors and sample_counts must be equal 1D arrays")
    if values.size == 0 or np.any(~np.isfinite(values)) or np.any(errors < 0.0):
        raise ValueError("batch estimates must be finite with nonnegative errors")
    if np.any(counts <= 0.0) or np.any(counts != np.floor(counts)):
        raise ValueError("sample counts must be positive integers")
    weights = counts / counts.sum()
    value = float(np.dot(weights, values))
    within = float(np.sqrt(np.dot(weights**2, errors**2)))
    between = None
    if values.size > 1:
        denominator = counts.sum() - float(np.dot(counts, counts) / counts.sum())
        variance = float(np.dot(counts, (values - value) ** 2) / denominator)
        effective_batches = float(counts.sum() ** 2 / np.dot(counts, counts))
        between = float(np.sqrt(variance / effective_batches))
    error = max(within, between if between is not None else 0.0)
    return CombinedEstimate(
        value=value,
        standard_error=error,
        within_batch_standard_error=within,
        between_batch_standard_error=between,
        total_trajectories=int(counts.sum()),
        batch_count=int(values.size),
    )


def assemble_eq40_rate(
    qtst_rates: np.ndarray,
    qtst_standard_errors: np.ndarray,
    channel_values: dict[tuple[int, int], float],
    channel_standard_errors: dict[tuple[int, int], float],
) -> tuple[float, float, list[tuple[int, int]]]:
    """Assemble Eq. (40), summing incoming states for each sampled PES.

    Keys are ``(surface_state, incoming_state)``.  Final electronic states are
    assumed already summed in each transmission batch.
    """

    rates = np.asarray(qtst_rates, dtype=float)
    rate_errors = np.asarray(qtst_standard_errors, dtype=float)
    if rates.ndim != 1 or rates.shape != rate_errors.shape:
        raise ValueError("QTST rates and errors must be equal 1D arrays")
    if np.any(rates < 0.0) or np.any(rate_errors < 0.0):
        raise ValueError("QTST rates and errors must be nonnegative")
    missing: list[tuple[int, int]] = []
    total = 0.0
    variance = 0.0
    for surface in range(rates.size):
        kappa_sum = 0.0
        kappa_variance = 0.0
        for incoming in range(rates.size):
            key = (surface, incoming)
            if key not in channel_values or key not in channel_standard_errors:
                missing.append(key)
                continue
            kappa_sum += float(channel_values[key])
            kappa_variance += float(channel_standard_errors[key]) ** 2
        total += rates[surface] * kappa_sum
        variance += (rate_errors[surface] * kappa_sum) ** 2
        variance += rates[surface] ** 2 * kappa_variance
    return float(total), float(np.sqrt(variance)), missing




