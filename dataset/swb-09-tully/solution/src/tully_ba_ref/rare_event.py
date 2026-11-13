"""Pure numerical kernels for Eqs. (50)-(56) of Shushkov et al. (2012)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def trial_hop_probability(projected_nac_rate: float, dt: float) -> float:
    r"""Equation (53): \(f_{ij}=\frac12[1-e^{-|\dot R d_{ij}|\Delta t}]\)."""

    if dt < 0.0:
        raise ValueError("dt must be nonnegative")
    rate = complex(projected_nac_rate)
    if abs(rate.imag) > 1.0e-12:
        raise ValueError("projected NAC rate must be real")
    return float(0.5 * (-np.expm1(-abs(rate.real) * dt)))


def trial_hop_probabilities(
    time_derivative_coupling: np.ndarray, active_state: int, dt: float
) -> np.ndarray:
    coupling = np.asarray(time_derivative_coupling)
    if coupling.ndim != 2 or coupling.shape[0] != coupling.shape[1]:
        raise ValueError("time_derivative_coupling must be a square matrix")
    probabilities = np.zeros(coupling.shape[0], dtype=float)
    for state in range(coupling.shape[0]):
        if state != active_state:
            probabilities[state] = trial_hop_probability(
                coupling[state, active_state], dt
            )
    if probabilities.sum() >= 1.0:
        raise ValueError("trial hop probability budget is not below one")
    return probabilities


def fssh_probabilities(
    rho: np.ndarray,
    time_derivative_coupling: np.ndarray,
    active_state: int,
    dt: float,
) -> np.ndarray:
    """Fewest-switches probabilities used as ``g_ij`` in Eqs. (55)-(56)."""

    rho = np.asarray(rho, dtype=complex)
    coupling = np.asarray(time_derivative_coupling, dtype=complex)
    if rho.shape != coupling.shape or rho.ndim != 2 or rho.shape[0] != rho.shape[1]:
        raise ValueError("rho and coupling must be equal square matrices")
    population = float(np.real(rho[active_state, active_state]))
    if population <= 0.0:
        raise ValueError("active-state population must be positive")
    bij = -2.0 * np.real(np.conj(rho) * coupling)
    probabilities = np.maximum(0.0, bij[:, active_state] * dt / population)
    probabilities[active_state] = 0.0
    if probabilities.sum() > 1.0 + 1.0e-12:
        raise ValueError("FSSH probability budget exceeds one; reduce dt")
    return probabilities


def backward_step_weight(
    fssh_probability: np.ndarray,
    trial_probability: np.ndarray,
    attempted_state: int | None,
) -> float:
    """Return one likelihood-ratio factor ``w_mu`` from Eq. (56)."""

    g = np.asarray(fssh_probability, dtype=float)
    f = np.asarray(trial_probability, dtype=float)
    if g.shape != f.shape or g.ndim != 1:
        raise ValueError("g and f must be equal 1D probability arrays")
    if np.any(g < 0.0) or np.any(f < 0.0) or g.sum() > 1.0 + 1.0e-12 or f.sum() >= 1.0:
        raise ValueError("invalid probability budget")
    if attempted_state is None:
        return float((1.0 - g.sum()) / (1.0 - f.sum()))
    if f[attempted_state] <= 0.0:
        raise ValueError("attempted event has zero trial probability")
    return float(g[attempted_state] / f[attempted_state])


def backward_path_weight(step_weights: np.ndarray) -> tuple[float, float]:
    """Equation (55), returned as both ``W`` and ``log(W)`` for diagnostics."""

    weights = np.asarray(step_weights, dtype=float)
    if weights.ndim != 1 or np.any(weights < 0.0):
        raise ValueError("step weights must be a nonnegative 1D array")
    if np.any(weights == 0.0):
        return 0.0, -np.inf
    log_weight = float(np.sum(np.log(weights)))
    if log_weight > np.log(np.finfo(float).max):
        return np.inf, log_weight
    if log_weight < np.log(np.finfo(float).tiny):
        return 0.0, log_weight
    return float(np.exp(log_weight)), log_weight


def reconstruct_dividing_surface_amplitudes(
    backward_propagator: np.ndarray,
    incoming_state: int,
    *,
    unitarity_tolerance: float = 1.0e-8,
) -> np.ndarray:
    r"""Recover \(c(0)\) from two independent backward solutions and Eq. (54).

    Columns of ``backward_propagator`` are the independent electronic
    solutions propagated from the dividing surface to the reactant asymptote.
    """

    propagator = np.asarray(backward_propagator, dtype=complex)
    if propagator.ndim != 2 or propagator.shape[0] != propagator.shape[1]:
        raise ValueError("backward_propagator must be square")
    identity = np.eye(propagator.shape[0], dtype=complex)
    defect = np.linalg.norm(propagator.conj().T @ propagator - identity, ord=2)
    if defect > unitarity_tolerance:
        raise ValueError(
            "backward_propagator is not unitary within the requested tolerance"
        )
    target = np.zeros(propagator.shape[0], dtype=complex)
    target[incoming_state] = 1.0
    amplitudes = np.linalg.solve(propagator, target)
    norm = np.linalg.norm(amplitudes)
    if norm <= 0.0:
        raise ValueError("reconstructed amplitudes have zero norm")
    return amplitudes / norm


def count_forward_crossings(
    centroid_positions: np.ndarray, dividing_surface: float, tol: float = 1.0e-12
) -> int:
    """Count reactant-to-product crossings, including departure from ``R=R‡``."""

    positions = np.asarray(centroid_positions, dtype=float)
    if positions.ndim != 1 or positions.size < 2:
        raise ValueError("at least two centroid positions are required")
    shifted = positions - float(dividing_surface)
    return int(np.count_nonzero((shifted[:-1] <= tol) & (shifted[1:] > tol)))


def effective_sample_size(weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=float)
    if weights.ndim != 1 or np.any(weights < 0.0):
        raise ValueError("weights must be a nonnegative 1D array")
    denominator = float(np.dot(weights, weights))
    if denominator == 0.0:
        return 0.0
    return float(weights.sum() ** 2 / denominator)


@dataclass(frozen=True)
class RatioEstimate:
    value: float
    standard_error: float
    block_values: np.ndarray
    effective_sample_size: float
    contribution_effective_sample_size: float


@dataclass(frozen=True)
class RatioTimeSeriesEstimate:
    """Block-resolved running estimate of the Appendix Eq. (50) ratio."""

    values: np.ndarray
    standard_errors: np.ndarray
    block_values: np.ndarray
    contribution_effective_sample_sizes: np.ndarray


def transmission_ratio_with_block_error(
    centroid_momenta: np.ndarray,
    mass: float,
    reactive_indicator: np.ndarray,
    backward_weights: np.ndarray,
    crossing_counts: np.ndarray,
    *,
    nblocks: int = 10,
) -> RatioEstimate:
    """Estimate Eq. (50) using positive outward flux and crossing correction."""

    pbar = np.asarray(centroid_momenta, dtype=float)
    reactive = np.asarray(reactive_indicator, dtype=float)
    weights = np.asarray(backward_weights, dtype=float)
    crossings = np.asarray(crossing_counts, dtype=int)
    if not (pbar.shape == reactive.shape == weights.shape == crossings.shape):
        raise ValueError("all trajectory arrays must have equal shape")
    if pbar.ndim != 1 or pbar.size < nblocks or nblocks < 2:
        raise ValueError("need at least one sample per block and at least two blocks")
    if (
        np.any(pbar <= 0.0)
        or np.any(weights < 0.0)
        or not np.all(np.isfinite(weights))
        or np.any(crossings < 0)
    ):
        raise ValueError("requires positive outward momenta and nonnegative weights/crossings")
    if np.any((reactive > 0.0) & (crossings < 1)):
        raise ValueError("every reactive trajectory must have a forward crossing")
    flux = pbar / float(mass)
    corrected = np.where(
        reactive > 0.0, reactive * weights / np.maximum(crossings, 1), 0.0
    )
    numerator = flux * corrected
    value = float(np.mean(numerator) / np.mean(flux))

    blocks = np.array_split(np.arange(pbar.size), nblocks)
    block_values = np.array(
        [np.mean(numerator[idx]) / np.mean(flux[idx]) for idx in blocks], dtype=float
    )
    standard_error = float(np.std(block_values, ddof=1) / np.sqrt(nblocks))
    return RatioEstimate(
        value=value,
        standard_error=standard_error,
        block_values=block_values,
        effective_sample_size=effective_sample_size(weights),
        contribution_effective_sample_size=effective_sample_size(numerator),
    )


def transmission_time_series_with_block_error(
    centroid_momenta: np.ndarray,
    mass: float,
    product_side: np.ndarray,
    backward_weights: np.ndarray,
    cumulative_forward_crossings: np.ndarray,
    *,
    nblocks: int = 10,
) -> RatioTimeSeriesEstimate:
    r"""Estimate a running flux-side ratio compatible with the Eq. (50) endpoint.

    ``product_side[a, t]`` is :math:`h(\bar R_a(t)-R^\ddagger)` and
    ``cumulative_forward_crossings[a, t]`` counts positive-direction departures
    through the dividing surface up to that time.  The inverse crossing count is
    the positive-outward-momentum correction described immediately after
    Eq. (56).  Consequently, the last time point is algebraically identical to
    :func:`transmission_ratio_with_block_error` when every trajectory has reached
    an asymptotic boundary.

    This is an operational time-resolved diagnostic for the paper's rare-event
    sampler.  It must not be confused with the separate RPMDrate ``+p/-p``
    signed-flux estimator.
    """

    pbar = np.asarray(centroid_momenta, dtype=float)
    side = np.asarray(product_side, dtype=bool)
    weights = np.asarray(backward_weights, dtype=float)
    crossings = np.asarray(cumulative_forward_crossings, dtype=int)
    if pbar.ndim != 1 or weights.shape != pbar.shape:
        raise ValueError("centroid momenta and weights must be equal 1D arrays")
    if side.ndim != 2 or crossings.shape != side.shape or side.shape[0] != pbar.size:
        raise ValueError("side and crossing histories must be equal trajectory-by-time arrays")
    if pbar.size < nblocks or nblocks < 2:
        raise ValueError("need at least one trajectory per block and at least two blocks")
    if (
        np.any(pbar <= 0.0)
        or np.any(weights < 0.0)
        or not np.all(np.isfinite(weights))
        or np.any(crossings < 0)
    ):
        raise ValueError("requires positive momenta and finite nonnegative weights/crossings")
    if np.any(side & (crossings < 1)):
        raise ValueError("a product-side sample must have at least one forward crossing")

    flux = pbar / float(mass)
    corrected = np.where(
        side,
        weights[:, None] / np.maximum(crossings, 1),
        0.0,
    )
    numerator = flux[:, None] * corrected
    denominator = float(np.mean(flux))
    values = np.mean(numerator, axis=0) / denominator

    blocks = np.array_split(np.arange(pbar.size), nblocks)
    block_values = np.asarray(
        [
            np.mean(numerator[idx], axis=0) / float(np.mean(flux[idx]))
            for idx in blocks
        ]
    )
    standard_errors = np.std(block_values, axis=0, ddof=1) / np.sqrt(nblocks)
    contribution_ess = np.asarray(
        [effective_sample_size(numerator[:, time]) for time in range(side.shape[1])]
    )
    return RatioTimeSeriesEstimate(
        values=np.asarray(values, dtype=float),
        standard_errors=np.asarray(standard_errors, dtype=float),
        block_values=np.asarray(block_values, dtype=float),
        contribution_effective_sample_sizes=contribution_ess,
    )




