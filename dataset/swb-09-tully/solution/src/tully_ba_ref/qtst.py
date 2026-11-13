"""Constrained-centroid path-integral mean-force estimator for MPES-QTST."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
try:
    from scipy.integrate import cumulative_simpson
except ImportError:  # SciPy < 1.12, used by the frozen jinsi environment
    from scipy.integrate import simpson

    def cumulative_simpson(values, *, x, initial=0.0):
        """Compatibility fallback; production on current SciPy uses the native routine."""

        values = np.asarray(values, dtype=float)
        grid = np.asarray(x, dtype=float)
        result = np.empty(values.size, dtype=float)
        result[0] = float(initial)
        for stop in range(2, values.size + 1):
            result[stop - 1] = float(initial) + float(
                simpson(values[:stop], x=grid[:stop])
            )
        return result

from .constants import KB_AU_PER_K as kB

from .trajectory import _spring_modes


@dataclass(frozen=True)
class MeanForceChain:
    centroid: float
    state: int
    mean_derivative: float
    standard_error: float
    effective_sample_size: float
    autocorrelation_time: float
    acceptance_rate: float
    proposal_scale: float
    nsamples: int
    series: np.ndarray


@dataclass(frozen=True)
class QTSTEstimate:
    grid: np.ndarray
    mean_force_by_seed: np.ndarray
    free_energy_by_seed: np.ndarray
    rates_by_seed: np.ndarray
    rate_mean: np.ndarray
    rate_standard_error: np.ndarray
    minimum_chain_ess: float
    acceptance_range: tuple[float, float]


def _autocorrelation_diagnostics(series: np.ndarray) -> tuple[float, float, float]:
    values = np.asarray(series, dtype=float)
    if values.ndim != 1 or values.size < 2:
        raise ValueError("series must contain at least two samples")
    if np.allclose(values, values[0], rtol=0.0, atol=1.0e-14):
        return 1.0, float(values.size), 0.0
    centered = values - np.mean(values)
    variance = float(np.dot(centered, centered) / values.size)
    nfft = 1 << (2 * values.size - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=nfft)
    autocovariance = np.fft.irfft(spectrum * spectrum.conj(), n=nfft)[: values.size]
    autocovariance /= np.arange(values.size, 0, -1)
    correlation = autocovariance / autocovariance[0]
    tau = 1.0
    for lag in range(1, values.size - 1, 2):
        pair = float(correlation[lag] + correlation[lag + 1])
        if pair <= 0.0:
            break
        tau += 2.0 * pair
    tau = min(max(tau, 1.0), float(values.size))
    ess = float(values.size / tau)
    standard_error = float(np.sqrt(variance / ess))
    return tau, ess, standard_error


def sample_constrained_mean_force(
    model,
    *,
    centroid: float,
    state: int,
    temperature: float,
    nbeads: int,
    nsteps: int,
    burnin_fraction: float,
    thin: int,
    seed: int,
    proposal_scale: float = 0.25,
    adapt_interval: int = 100,
    target_acceptance: float = 0.3,
) -> MeanForceChain:
    r"""Sample a fixed-centroid path with a free-ring-polymer pCN proposal.

    The pCN move is reversible with respect to the exact constrained free-ring
    Gaussian, so the Metropolis ratio contains only
    \(\exp[-\beta_P\Delta\sum_k V_s(q_k)]\).  This avoids the severe random-walk
    scaling of one-bead moves at 128 beads.
    """

    if (
        temperature <= 0.0
        or nbeads < 2
        or nsteps < 2
        or not 0.0 <= burnin_fraction < 1.0
        or thin < 1
        or not 0.0 < proposal_scale < 1.0
    ):
        raise ValueError("invalid constrained-chain parameters")
    burnin = int(nsteps * burnin_fraction)
    if nsteps - burnin < 2 * thin:
        raise ValueError("too few production samples")
    beta = 1.0 / (kB * float(temperature))
    beta_p = beta / nbeads
    mass = float(model.mass)
    eigenvalues, modes = _spring_modes(nbeads)
    mode_std = np.zeros(nbeads)
    nonzero = eigenvalues > 1.0e-14
    mode_std[nonzero] = np.sqrt(beta_p / (mass * eigenvalues[nonzero]))
    rng = np.random.default_rng(seed)
    path = np.full(nbeads, float(centroid))

    def potential_sum(candidate):
        return float(np.sum(model.get_energies(candidate)[:, state]))

    potential = potential_sum(path)
    logit = float(np.log(proposal_scale / (1.0 - proposal_scale)))
    accepted_total = 0
    accepted_window = 0
    proposed_window = 0
    derivatives = []
    for step in range(nsteps):
        scale = 1.0 / (1.0 + np.exp(-logit))
        gaussian_modes = mode_std * rng.normal(size=nbeads)
        gaussian_modes[~nonzero] = 0.0
        gaussian_internal = modes @ gaussian_modes
        current_internal = path - float(centroid)
        candidate = (
            float(centroid)
            + np.sqrt(1.0 - scale**2) * current_internal
            + scale * gaussian_internal
        )
        candidate += float(centroid) - np.mean(candidate)
        candidate_potential = potential_sum(candidate)
        log_ratio = -beta_p * (candidate_potential - potential)
        accepted = log_ratio >= 0.0 or np.log(rng.random()) < log_ratio
        if accepted:
            path = candidate
            potential = candidate_potential
            accepted_total += 1
        if step < burnin:
            proposed_window += 1
            accepted_window += int(accepted)
            if proposed_window == adapt_interval:
                observed = accepted_window / proposed_window
                logit += 0.5 * (observed - target_acceptance)
                logit = float(np.clip(logit, -4.595, 4.595))
                accepted_window = 0
                proposed_window = 0
        elif (step - burnin) % thin == 0:
            force = model.evaluate(
                path, need_nac=False, need_wavefun=False
            ).force_diag[:, state]
            derivatives.append(float(-np.mean(force)))
    series = np.asarray(derivatives)
    tau, ess, standard_error = _autocorrelation_diagnostics(series)
    return MeanForceChain(
        centroid=float(centroid),
        state=int(state),
        mean_derivative=float(np.mean(series)),
        standard_error=standard_error,
        effective_sample_size=ess,
        autocorrelation_time=tau,
        acceptance_rate=float(accepted_total / nsteps),
        proposal_scale=float(1.0 / (1.0 + np.exp(-logit))),
        nsamples=int(series.size),
        series=series,
    )


def integrate_mean_force(grid: np.ndarray, derivative: np.ndarray) -> np.ndarray:
    """Integrate a physical centroid PMF derivative with cumulative Simpson quadrature."""

    grid = np.asarray(grid, dtype=float)
    derivative = np.asarray(derivative, dtype=float)
    if grid.ndim != 1 or derivative.shape != grid.shape or grid.size < 2:
        raise ValueError("grid and derivative must be equal 1D arrays")
    if np.any(np.diff(grid) <= 0.0):
        raise ValueError("grid must be strictly increasing")
    return cumulative_simpson(derivative, x=grid, initial=0.0)


def estimate_mpes_qtst(
    model,
    *,
    temperature: float,
    nbeads: int,
    reactant_surface: float,
    dividing_surface: float,
    nwindows: int,
    nsteps: int,
    seeds: tuple[int, ...],
    burnin_fraction: float = 0.25,
    thin: int = 5,
    proposal_scale: float = 0.25,
) -> tuple[QTSTEstimate, list[MeanForceChain]]:
    """Estimate the two state-selected Eq. (42) MPES-QTST rates."""

    if nwindows < 2 or len(seeds) < 2:
        raise ValueError("need at least two windows and two independent seeds")
    grid = np.linspace(reactant_surface, dividing_surface, nwindows)
    chains: list[MeanForceChain] = []
    derivatives = np.empty((len(seeds), model.nstate, nwindows))
    acceptances = []
    ess_values = []
    for seed_index, seed in enumerate(seeds):
        for state in range(model.nstate):
            for window_index, centroid in enumerate(grid):
                chain = sample_constrained_mean_force(
                    model,
                    centroid=float(centroid),
                    state=state,
                    temperature=temperature,
                    nbeads=nbeads,
                    nsteps=nsteps,
                    burnin_fraction=burnin_fraction,
                    thin=thin,
                    seed=int(seed + 10000 * state + window_index),
                    proposal_scale=proposal_scale,
                )
                chains.append(chain)
                derivatives[seed_index, state, window_index] = chain.mean_derivative
                acceptances.append(chain.acceptance_rate)
                ess_values.append(chain.effective_sample_size)
    free_energy = np.empty_like(derivatives)
    for seed_index in range(len(seeds)):
        for state in range(model.nstate):
            free_energy[seed_index, state] = integrate_mean_force(
                grid, derivatives[seed_index, state]
            )
    beta = 1.0 / (kB * float(temperature))
    free_flux = 0.5 / np.sqrt(2.0 * np.pi * float(model.mass) * beta)
    rates = free_flux * np.exp(-beta * free_energy[:, :, -1])
    return (
        QTSTEstimate(
            grid=grid,
            mean_force_by_seed=derivatives,
            free_energy_by_seed=free_energy,
            rates_by_seed=rates,
            rate_mean=np.mean(rates, axis=0),
            rate_standard_error=np.std(rates, axis=0, ddof=1)
            / np.sqrt(len(seeds)),
            minimum_chain_ess=float(np.min(ess_values)),
            acceptance_range=(float(np.min(acceptances)), float(np.max(acceptances))),
        ),
        chains,
    )




