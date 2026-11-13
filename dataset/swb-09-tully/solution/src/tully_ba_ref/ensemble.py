"""Small-to-large statistical batches for the Eq. (50) transmission coefficient."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bidirectional import run_bidirectional_trajectory
from .rare_event import (
    RatioEstimate,
    RatioTimeSeriesEstimate,
    transmission_ratio_with_block_error,
    transmission_time_series_with_block_error,
)
from .sampling import sample_positive_flux_momenta
from .trajectory import TrajectoryState


@dataclass(frozen=True)
class TransmissionBatch:
    total: RatioEstimate
    by_final_state: tuple[RatioEstimate, ...]
    centroid_momenta: np.ndarray
    backward_weights: np.ndarray
    reactive: np.ndarray
    final_states: np.ndarray
    crossing_counts: np.ndarray
    backward_attempts: np.ndarray
    backward_accepted: np.ndarray
    forward_attempts: np.ndarray
    forward_accepted: np.ndarray


@dataclass(frozen=True)
class FluxSideBatch:
    """Time-resolved Appendix rare-event diagnostic on a common time grid."""

    time_steps: np.ndarray
    total: RatioTimeSeriesEstimate
    by_active_state: tuple[RatioTimeSeriesEstimate, ...]
    endpoint: TransmissionBatch


def run_transmission_batch(
    model,
    configurations: np.ndarray,
    *,
    surface_state: int,
    incoming_state: int,
    temperature: float,
    dt: float,
    approximation: str,
    dividing_surface: float,
    reactant_boundary: float,
    product_boundary: float,
    max_steps: int,
    seed: int,
    nblocks: int = 10,
    event_convention: str = "stored_backward",
    frustrated_hop_policy: str = "keep",
    electronic_integrator: str = "midpoint_exponential",
) -> TransmissionBatch:
    """Run a reproducible batch without retaining bead/electronic time histories."""

    configurations = np.asarray(configurations, dtype=float)
    if configurations.ndim != 2 or configurations.shape[0] < nblocks:
        raise ValueError("configurations must contain at least one trajectory per block")
    ntraj, nbeads = configurations.shape
    pbar = np.empty(ntraj)
    weights = np.empty(ntraj)
    reactive = np.empty(ntraj, dtype=bool)
    final_states = np.empty(ntraj, dtype=int)
    crossings = np.empty(ntraj, dtype=int)
    backward_attempts = np.empty(ntraj, dtype=int)
    backward_accepted = np.empty(ntraj, dtype=int)
    forward_attempts = np.empty(ntraj, dtype=int)
    forward_accepted = np.empty(ntraj, dtype=int)
    seed_sequence = np.random.SeedSequence(seed).spawn(ntraj)
    for trajectory, (positions, child_seed) in enumerate(zip(configurations, seed_sequence)):
        momentum_rng, backward_rng, forward_rng = [
            np.random.default_rng(item) for item in child_seed.spawn(3)
        ]
        momenta = sample_positive_flux_momenta(
            mass=float(model.mass),
            temperature=temperature,
            nbeads=nbeads,
            rng=momentum_rng,
        )
        amplitudes = np.zeros(model.nstate, dtype=complex)
        amplitudes[surface_state] = 1.0
        initial = TrajectoryState(
            positions=positions.copy(),
            momenta=momenta,
            active_state=surface_state,
            amplitudes=amplitudes,
        )
        try:
            result = run_bidirectional_trajectory(
                model,
                initial,
                incoming_state=incoming_state,
                temperature=temperature,
                dt=dt,
                approximation=approximation,
                dividing_surface=dividing_surface,
                reactant_boundary=reactant_boundary,
                product_boundary=product_boundary,
                max_steps=max_steps,
                backward_rng=backward_rng,
                forward_rng=forward_rng,
                nuclear_substeps=1,
                compute_step_energies=False,
                event_convention=event_convention,
                frustrated_hop_policy=frustrated_hop_policy,
                electronic_integrator=electronic_integrator,
            )
        except RuntimeError as error:
            raise RuntimeError(f"trajectory {trajectory} did not complete: {error}") from error
        pbar[trajectory] = float(np.mean(momenta))
        weights[trajectory] = (
            result.reweighting.weight if result.reweighting is not None else 0.0
        )
        reactive[trajectory] = result.reactive
        final_states[trajectory] = result.final_state
        crossings[trajectory] = result.forward_crossings
        backward_attempts[trajectory] = sum(
            item is not None for item in result.backward.attempted_states
        )
        backward_accepted[trajectory] = int(np.sum(result.backward.hop_acceptances))
        forward_attempts[trajectory] = (
            sum(item is not None for item in result.forward.attempted_states)
            if result.forward is not None
            else 0
        )
        forward_accepted[trajectory] = (
            int(np.sum(result.forward.hop_acceptances))
            if result.forward is not None
            else 0
        )
    total = transmission_ratio_with_block_error(
        pbar,
        float(model.mass),
        reactive,
        weights,
        crossings,
        nblocks=nblocks,
    )
    by_final_state = tuple(
        transmission_ratio_with_block_error(
            pbar,
            float(model.mass),
            reactive & (final_states == state),
            weights,
            crossings,
            nblocks=nblocks,
        )
        for state in range(model.nstate)
    )
    return TransmissionBatch(
        total=total,
        by_final_state=by_final_state,
        centroid_momenta=pbar,
        backward_weights=weights,
        reactive=reactive,
        final_states=final_states,
        crossing_counts=crossings,
        backward_attempts=backward_attempts,
        backward_accepted=backward_accepted,
        forward_attempts=forward_attempts,
        forward_accepted=forward_accepted,
    )


def run_flux_side_batch(
    model,
    configurations: np.ndarray,
    *,
    surface_state: int,
    incoming_state: int,
    temperature: float,
    dt: float,
    approximation: str,
    dividing_surface: float,
    reactant_boundary: float,
    product_boundary: float,
    max_steps: int,
    seed: int,
    sample_stride: int = 100,
    nblocks: int = 10,
    event_convention: str = "stored_backward",
    frustrated_hop_policy: str = "keep",
    electronic_integrator: str = "midpoint_exponential",
) -> FluxSideBatch:
    """Retain a compact running flux-side history while reproducing Eq. (50).

    Trajectories are propagated exactly as in :func:`run_transmission_batch`.
    Once an asymptotic boundary is reached, its side and electronic state are
    held fixed on the remaining common time grid.  This is valid for the present
    scattering boundaries, where the NAC has vanished and the outgoing path
    cannot turn back.
    """

    configurations = np.asarray(configurations, dtype=float)
    if configurations.ndim != 2 or configurations.shape[0] < nblocks:
        raise ValueError("configurations must contain at least one trajectory per block")
    if sample_stride < 1:
        raise ValueError("sample_stride must be positive")
    ntraj, nbeads = configurations.shape
    sample_steps = np.unique(
        np.concatenate(
            [np.arange(1, max_steps + 1, sample_stride, dtype=int), [max_steps]]
        )
    )
    ntime = sample_steps.size
    product_side = np.zeros((ntraj, ntime), dtype=bool)
    active_states = np.full((ntraj, ntime), -1, dtype=int)
    cumulative_crossings = np.zeros((ntraj, ntime), dtype=int)
    maximum_forward_steps = 0

    pbar = np.empty(ntraj)
    weights = np.empty(ntraj)
    reactive = np.empty(ntraj, dtype=bool)
    final_states = np.empty(ntraj, dtype=int)
    crossings = np.empty(ntraj, dtype=int)
    backward_attempts = np.empty(ntraj, dtype=int)
    backward_accepted = np.empty(ntraj, dtype=int)
    forward_attempts = np.empty(ntraj, dtype=int)
    forward_accepted = np.empty(ntraj, dtype=int)
    seed_sequence = np.random.SeedSequence(seed).spawn(ntraj)

    for trajectory, (positions, child_seed) in enumerate(zip(configurations, seed_sequence)):
        momentum_rng, backward_rng, forward_rng = [
            np.random.default_rng(item) for item in child_seed.spawn(3)
        ]
        momenta = sample_positive_flux_momenta(
            mass=float(model.mass),
            temperature=temperature,
            nbeads=nbeads,
            rng=momentum_rng,
        )
        amplitudes = np.zeros(model.nstate, dtype=complex)
        amplitudes[surface_state] = 1.0
        initial = TrajectoryState(
            positions=positions.copy(),
            momenta=momenta,
            active_state=surface_state,
            amplitudes=amplitudes,
        )
        try:
            result = run_bidirectional_trajectory(
                model,
                initial,
                incoming_state=incoming_state,
                temperature=temperature,
                dt=dt,
                approximation=approximation,
                dividing_surface=dividing_surface,
                reactant_boundary=reactant_boundary,
                product_boundary=product_boundary,
                max_steps=max_steps,
                backward_rng=backward_rng,
                forward_rng=forward_rng,
                nuclear_substeps=1,
                compute_step_energies=False,
                event_convention=event_convention,
                frustrated_hop_policy=frustrated_hop_policy,
                electronic_integrator=electronic_integrator,
            )
        except RuntimeError as error:
            raise RuntimeError(f"trajectory {trajectory} did not complete: {error}") from error

        pbar[trajectory] = float(np.mean(momenta))
        weights[trajectory] = (
            result.reweighting.weight if result.reweighting is not None else 0.0
        )
        reactive[trajectory] = result.reactive
        final_states[trajectory] = result.final_state
        crossings[trajectory] = result.forward_crossings
        backward_attempts[trajectory] = sum(
            item is not None for item in result.backward.attempted_states
        )
        backward_accepted[trajectory] = int(np.sum(result.backward.hop_acceptances))
        forward_attempts[trajectory] = (
            sum(item is not None for item in result.forward.attempted_states)
            if result.forward is not None
            else 0
        )
        forward_accepted[trajectory] = (
            int(np.sum(result.forward.hop_acceptances))
            if result.forward is not None
            else 0
        )

        if result.forward is None:
            continue
        forward_positions = np.asarray(result.forward.centroid_positions, dtype=float)
        maximum_forward_steps = max(maximum_forward_steps, forward_positions.size - 1)
        forward_states = np.asarray(result.forward.active_states, dtype=int)
        shifted = forward_positions - float(dividing_surface)
        departures = (shifted[:-1] <= 1.0e-12) & (shifted[1:] > 1.0e-12)
        crossing_history = np.concatenate([[0], np.cumsum(departures, dtype=int)])
        indices = np.minimum(sample_steps, forward_positions.size - 1)
        sampled_positions = forward_positions[indices]
        product_side[trajectory] = sampled_positions > float(dividing_surface) + 1.0e-12
        active_states[trajectory] = forward_states[indices]
        cumulative_crossings[trajectory] = crossing_history[indices]

    # The propagation cap is deliberately generous.  Keep only the common grid
    # through the first sampled point at or after the longest realized forward
    # path; all shorter scattering paths have already been extended by their
    # asymptotic side/state.  This avoids thousands of redundant plateau points.
    keep_points = max(
        1,
        min(
            ntime,
            int(np.searchsorted(sample_steps, maximum_forward_steps, side="left")) + 1,
        ),
    )
    sample_steps = sample_steps[:keep_points]
    product_side = product_side[:, :keep_points]
    active_states = active_states[:, :keep_points]
    cumulative_crossings = cumulative_crossings[:, :keep_points]

    endpoint_total = transmission_ratio_with_block_error(
        pbar, float(model.mass), reactive, weights, crossings, nblocks=nblocks
    )
    endpoint_by_state = tuple(
        transmission_ratio_with_block_error(
            pbar,
            float(model.mass),
            reactive & (final_states == state),
            weights,
            crossings,
            nblocks=nblocks,
        )
        for state in range(model.nstate)
    )
    endpoint = TransmissionBatch(
        total=endpoint_total,
        by_final_state=endpoint_by_state,
        centroid_momenta=pbar,
        backward_weights=weights,
        reactive=reactive,
        final_states=final_states,
        crossing_counts=crossings,
        backward_attempts=backward_attempts,
        backward_accepted=backward_accepted,
        forward_attempts=forward_attempts,
        forward_accepted=forward_accepted,
    )
    total = transmission_time_series_with_block_error(
        pbar,
        float(model.mass),
        product_side,
        weights,
        cumulative_crossings,
        nblocks=nblocks,
    )
    by_active_state = tuple(
        transmission_time_series_with_block_error(
            pbar,
            float(model.mass),
            product_side & (active_states == state),
            weights,
            cumulative_crossings,
            nblocks=nblocks,
        )
        for state in range(model.nstate)
    )
    if not np.isclose(total.values[-1], endpoint.total.value, rtol=0.0, atol=1.0e-14):
        raise RuntimeError("running flux-side endpoint does not match Eq. (50) estimator")
    return FluxSideBatch(
        time_steps=sample_steps,
        total=total,
        by_active_state=by_active_state,
        endpoint=endpoint,
    )




