"""Backward trial-path replay and the Eqs. (54)-(56) likelihood ratio."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .rare_event import (
    backward_path_weight,
    backward_step_weight,
    fssh_probabilities,
    reconstruct_dividing_surface_amplitudes,
    trial_hop_probabilities,
)


@dataclass(frozen=True)
class BackwardRecord:
    propagators: np.ndarray
    time_derivative_couplings: np.ndarray
    active_states: np.ndarray
    attempted_states: tuple[int | None, ...]
    dt: float
    active_states_after: np.ndarray | None = None
    hop_acceptances: np.ndarray | None = None


@dataclass(frozen=True)
class BackwardReweighting:
    dividing_surface_amplitudes: np.ndarray
    amplitude_history: np.ndarray
    step_weights: np.ndarray
    weight: float
    log_weight: float


def reconstruct_backward_amplitude_history(
    propagators: np.ndarray, incoming_state: int
) -> np.ndarray:
    """Apply Eq. (54) to a stored sequence running from 0 to negative time."""

    steps = np.asarray(propagators, dtype=complex)
    if steps.ndim != 3 or steps.shape[1] != steps.shape[2]:
        raise ValueError("propagators must have shape (nstep, nstate, nstate)")
    total = np.eye(steps.shape[1], dtype=complex)
    for propagator in steps:
        total = propagator @ total
    c0 = reconstruct_dividing_surface_amplitudes(total, incoming_state)
    history = [c0]
    amplitude = c0
    for propagator in steps:
        amplitude = propagator @ amplitude
        history.append(amplitude)
    return np.asarray(history)


def reweight_backward_record(
    record: BackwardRecord,
    incoming_state: int,
    *,
    event_convention: str = "stored_backward",
) -> BackwardReweighting:
    """Reconstruct amplitudes and evaluate every Eq. (56) path factor."""

    couplings = np.asarray(record.time_derivative_couplings, dtype=complex)
    active = np.asarray(record.active_states, dtype=int)
    nsteps = np.asarray(record.propagators).shape[0]
    if couplings.shape[0] != nsteps or active.shape != (nsteps,):
        raise ValueError("record arrays disagree on the number of steps")
    if len(record.attempted_states) != nsteps or record.dt <= 0.0:
        raise ValueError("attempt history or dt is invalid")
    if event_convention not in {
        "stored_backward",
        "backward_replay",
        "backward_replay_newer",
        "physical_forward_reverse",
        "physical_forward_bijective",
        "physical_forward_bijective_negative",
    }:
        raise ValueError(
            "event_convention must be 'stored_backward', "
            "'physical_forward_reverse', or 'physical_forward_bijective'"
        )
    amplitudes = reconstruct_backward_amplitude_history(
        record.propagators, incoming_state
    )
    active_after = (
        None
        if record.active_states_after is None
        else np.asarray(record.active_states_after, dtype=int)
    )
    acceptances = (
        None
        if record.hop_acceptances is None
        else np.asarray(record.hop_acceptances, dtype=bool)
    )
    if event_convention in {
        "physical_forward_reverse",
        "physical_forward_bijective",
        "physical_forward_bijective_negative",
    }:
        if active_after is None or acceptances is None:
            raise ValueError(
                "physical-forward convention requires active_states_after and "
                "hop_acceptances"
            )
        if active_after.shape != (nsteps,) or acceptances.shape != (nsteps,):
            raise ValueError("physical reverse event arrays have invalid shape")
    factors = np.empty(nsteps, dtype=float)
    for step in range(nsteps):
        f = trial_hop_probabilities(couplings[step], active[step], record.dt)
        attempted_state = record.attempted_states[step]
        if event_convention == "stored_backward":
            # Historical implementation: evaluate the target probability in
            # the auxiliary backward direction using the newer-end active
            # state and the sign-reversed NAC velocity.
            rho = np.outer(
                amplitudes[step + 1], amplitudes[step + 1].conj()
            )
            g = fssh_probabilities(
                rho, -couplings[step], active[step], record.dt
            )
        elif event_convention in {"backward_replay", "backward_replay_newer"}:
            # Directly replay the event labels of the auxiliary backward path
            # with the physical-sign coupling.  The two diagnostic variants
            # isolate whether the hop test is attached to the older or newer
            # electronic endpoint.
            amplitude_index = (
                step + 1 if event_convention == "backward_replay" else step
            )
            rho = np.outer(
                amplitudes[amplitude_index], amplitudes[amplitude_index].conj()
            )
            g = fssh_probabilities(
                rho, couplings[step], active[step], record.dt
            )
        else:
            # A stored interval runs newer -> older.  Its physical forward
            # reverse starts from the older endpoint, with the physical NAC
            # velocity sign, propagates the amplitudes to the newer endpoint,
            # and only then applies the FSSH hop test.  An accepted backward
            # i->j proposal corresponds to a physical-forward j->i event.
            rho = np.outer(amplitudes[step], amplitudes[step].conj())
            physical_active = int(active_after[step])
            target_coupling = (
                -couplings[step]
                if event_convention == "physical_forward_bijective_negative"
                else couplings[step]
            )
            g = fssh_probabilities(rho, target_coupling, physical_active, record.dt)
            if (
                event_convention == "physical_forward_reverse"
                and attempted_state is not None
                and bool(acceptances[step])
            ):
                attempted_state = int(active[step])
        if (
            event_convention
            in {"physical_forward_bijective", "physical_forward_bijective_negative"}
            and attempted_state is not None
        ):
            proposal_target = int(attempted_state)
            if f[proposal_target] <= 0.0:
                raise ValueError("attempted proposal event has zero trial probability")
            physical_target = (
                int(active[step]) if bool(acceptances[step]) else proposal_target
            )
            factors[step] = float(g[physical_target] / f[proposal_target])
        else:
            factors[step] = backward_step_weight(g, f, attempted_state)
    weight, log_weight = backward_path_weight(factors)
    return BackwardReweighting(
        dividing_surface_amplitudes=amplitudes[0],
        amplitude_history=amplitudes,
        step_weights=factors,
        weight=weight,
        log_weight=log_weight,
    )




