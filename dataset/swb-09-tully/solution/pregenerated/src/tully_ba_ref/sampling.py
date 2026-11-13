"""Compact constrained configurations and positive-flux momenta for Eq. (50)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import KB_AU_PER_K as kB

from .trajectory import _spring_modes


@dataclass(frozen=True)
class ConstrainedConfigurations:
    positions: np.ndarray
    acceptance_rate: float
    final_proposal_scale: float
    burnin_steps: int
    steps_between: int


def sample_constrained_configurations(
    model,
    *,
    centroid: float,
    state: int,
    temperature: float,
    nbeads: int,
    nconfigurations: int,
    burnin_steps: int,
    steps_between: int,
    seed: int,
    proposal_scale: float = 0.5,
) -> ConstrainedConfigurations:
    """Draw fixed-centroid full-potential paths with a free-reference pCN kernel."""

    if nconfigurations < 1 or burnin_steps < 1 or steps_between < 1:
        raise ValueError("configuration counts and sampling steps must be positive")
    beta_p = 1.0 / (kB * temperature * nbeads)
    eigenvalues, modes = _spring_modes(nbeads)
    nonzero = eigenvalues > 1.0e-14
    mode_std = np.zeros(nbeads)
    mode_std[nonzero] = np.sqrt(beta_p / (float(model.mass) * eigenvalues[nonzero]))
    rng = np.random.default_rng(seed)
    path = np.full(nbeads, float(centroid))
    potential = float(np.sum(model.get_energies(path)[:, state]))
    logit = float(np.log(proposal_scale / (1.0 - proposal_scale)))
    accepted = 0
    total_steps = burnin_steps + nconfigurations * steps_between
    configurations = []
    for step in range(total_steps):
        scale = float(1.0 / (1.0 + np.exp(-logit)))
        modes_sample = mode_std * rng.normal(size=nbeads)
        modes_sample[~nonzero] = 0.0
        reference = modes @ modes_sample
        candidate = (
            float(centroid)
            + np.sqrt(1.0 - scale**2) * (path - float(centroid))
            + scale * reference
        )
        candidate += float(centroid) - np.mean(candidate)
        candidate_potential = float(np.sum(model.get_energies(candidate)[:, state]))
        log_ratio = -beta_p * (candidate_potential - potential)
        move_accepted = log_ratio >= 0.0 or np.log(rng.random()) < log_ratio
        if move_accepted:
            path = candidate
            potential = candidate_potential
            accepted += 1
        if step < burnin_steps and (step + 1) % 100 == 0:
            # Cumulative Robbins-style adaptation is used only before samples.
            observed = accepted / (step + 1)
            logit += 0.25 * (observed - 0.3) / np.sqrt((step + 1) / 100)
            logit = float(np.clip(logit, -4.595, 4.595))
        if step >= burnin_steps and (step - burnin_steps + 1) % steps_between == 0:
            configurations.append(path.copy())
    positions = np.asarray(configurations)
    if not np.allclose(np.mean(positions, axis=1), centroid, atol=1.0e-13):
        raise RuntimeError("constrained sampler lost the requested centroid")
    return ConstrainedConfigurations(
        positions=positions,
        acceptance_rate=float(accepted / total_steps),
        final_proposal_scale=float(1.0 / (1.0 + np.exp(-logit))),
        burnin_steps=burnin_steps,
        steps_between=steps_between,
    )


def sample_positive_flux_momenta(
    *, mass: float, temperature: float, nbeads: int, rng: np.random.Generator
) -> np.ndarray:
    """Draw canonical bead momenta conditioned on positive centroid momentum."""

    beta = 1.0 / (kB * temperature)
    momenta = rng.normal(0.0, np.sqrt(mass * nbeads / beta), size=nbeads)
    centroid = float(np.mean(momenta))
    if centroid < 0.0:
        momenta -= 2.0 * centroid
    return momenta




