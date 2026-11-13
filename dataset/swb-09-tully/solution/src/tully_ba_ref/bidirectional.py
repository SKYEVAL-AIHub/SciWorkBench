"""End-to-end backward/forward rare-event trajectory orchestration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .backward import BackwardRecord, BackwardReweighting, reweight_backward_record
from .rare_event import count_forward_crossings
from .trajectory import TrajectoryState, propagate_rpsh_step


@dataclass(frozen=True)
class PathLeg:
    centroid_positions: np.ndarray
    active_states: np.ndarray
    attempted_states: tuple[int | None, ...]
    hop_acceptances: np.ndarray
    completed: bool
    terminal_side: str
    final_state: TrajectoryState


@dataclass(frozen=True)
class BidirectionalResult:
    backward: PathLeg
    forward: PathLeg | None
    reweighting: BackwardReweighting | None
    incoming_state: int
    final_state: int
    reactive: bool
    forward_crossings: int
    corrected_weight: float


def run_backward_trial_path(
    model,
    dividing_surface_state: TrajectoryState,
    *,
    incoming_state: int,
    temperature: float,
    dt: float,
    approximation: str,
    reactant_boundary: float,
    product_boundary: float,
    max_steps: int,
    rng: np.random.Generator,
    nuclear_substeps: int = 1,
    compute_step_energies: bool = True,
    event_convention: str = "stored_backward",
    frustrated_hop_policy: str = "keep",
    electronic_integrator: str = "midpoint_exponential",
) -> tuple[PathLeg, BackwardReweighting]:
    """Run Eq. (53) trial hops backward and then apply Eqs. (54)-(56)."""

    if dt <= 0.0 or max_steps < 1:
        raise ValueError("dt and max_steps must be positive")
    state = dividing_surface_state
    centroids = [float(np.mean(state.positions))]
    active_before: list[int] = []
    active_history = [state.active_state]
    attempted: list[int | None] = []
    accepted: list[bool] = []
    propagators = []
    couplings = []
    completed = (
        centroids[-1] <= reactant_boundary or centroids[-1] >= product_boundary
    )
    while not completed and len(propagators) < max_steps:
        active_before.append(state.active_state)
        result = propagate_rpsh_step(
            model,
            state,
            temperature=temperature,
            dt=-dt,
            approximation=approximation,
            random_number=float(rng.random()),
            trial_hop=True,
            frustrated_hop_policy=frustrated_hop_policy,
            electronic_integrator=electronic_integrator,
            nuclear_substeps=nuclear_substeps,
            compute_energy=compute_step_energies,
        )
        propagators.append(result.propagator)
        couplings.append(result.time_derivative_coupling)
        attempted.append(result.attempted_state)
        accepted.append(result.hop_accepted)
        state = result.state
        centroids.append(float(np.mean(state.positions)))
        active_history.append(state.active_state)
        completed = (
            centroids[-1] <= reactant_boundary or centroids[-1] >= product_boundary
        )
    if not completed:
        raise RuntimeError(
            "backward path reached neither asymptotic boundary; "
            f"steps={len(propagators)}, final_centroid={centroids[-1]:.8g}, "
            f"range=[{min(centroids):.8g}, {max(centroids):.8g}]"
        )
    if not propagators:
        raise ValueError("backward path must contain at least one step")
    terminal_side = "reactant" if centroids[-1] <= reactant_boundary else "product"
    reweighting = None
    if terminal_side == "reactant" and state.active_state == incoming_state:
        record = BackwardRecord(
            propagators=np.asarray(propagators),
            time_derivative_couplings=np.asarray(couplings),
            active_states=np.asarray(active_before),
            attempted_states=tuple(attempted),
            dt=dt,
            active_states_after=np.asarray(active_history[1:]),
            hop_acceptances=np.asarray(accepted, dtype=bool),
        )
        reweighting = reweight_backward_record(
            record, incoming_state, event_convention=event_convention
        )
    leg = PathLeg(
        centroid_positions=np.asarray(centroids),
        active_states=np.asarray(active_history),
        attempted_states=tuple(attempted),
        hop_acceptances=np.asarray(accepted, dtype=bool),
        completed=completed,
        terminal_side=terminal_side,
        final_state=state,
    )
    return leg, reweighting


def run_forward_rpsh_path(
    model,
    dividing_surface_state: TrajectoryState,
    *,
    temperature: float,
    dt: float,
    approximation: str,
    reactant_boundary: float,
    product_boundary: float,
    max_steps: int,
    rng: np.random.Generator,
    nuclear_substeps: int = 1,
    compute_step_energies: bool = True,
    frustrated_hop_policy: str = "keep",
    electronic_integrator: str = "midpoint_exponential",
) -> PathLeg:
    """Run ordinary FSSH forward until either asymptotic boundary is reached."""

    if dt <= 0.0 or max_steps < 1 or reactant_boundary >= product_boundary:
        raise ValueError("invalid time step, step limit, or boundaries")
    state = dividing_surface_state
    centroids = [float(np.mean(state.positions))]
    active_history = [state.active_state]
    attempted: list[int | None] = []
    accepted: list[bool] = []
    completed = not (reactant_boundary < centroids[-1] < product_boundary)
    while not completed and len(attempted) < max_steps:
        result = propagate_rpsh_step(
            model,
            state,
            temperature=temperature,
            dt=dt,
            approximation=approximation,
            random_number=float(rng.random()),
            trial_hop=False,
            frustrated_hop_policy=frustrated_hop_policy,
            electronic_integrator=electronic_integrator,
            nuclear_substeps=nuclear_substeps,
            compute_energy=compute_step_energies,
        )
        attempted.append(result.attempted_state)
        accepted.append(result.hop_accepted)
        state = result.state
        centroids.append(float(np.mean(state.positions)))
        active_history.append(state.active_state)
        completed = not (reactant_boundary < centroids[-1] < product_boundary)
    return PathLeg(
        centroid_positions=np.asarray(centroids),
        active_states=np.asarray(active_history),
        attempted_states=tuple(attempted),
        hop_acceptances=np.asarray(accepted, dtype=bool),
        completed=completed,
        terminal_side=(
            "reactant"
            if centroids[-1] <= reactant_boundary
            else "product"
            if centroids[-1] >= product_boundary
            else "max_steps"
        ),
        final_state=state,
    )


def run_bidirectional_trajectory(
    model,
    dividing_surface_state: TrajectoryState,
    *,
    incoming_state: int,
    temperature: float,
    dt: float,
    approximation: str,
    dividing_surface: float,
    reactant_boundary: float,
    product_boundary: float,
    max_steps: int,
    backward_rng: np.random.Generator,
    forward_rng: np.random.Generator,
    nuclear_substeps: int = 1,
    compute_step_energies: bool = True,
    event_convention: str = "stored_backward",
    frustrated_hop_policy: str = "keep",
    electronic_integrator: str = "midpoint_exponential",
) -> BidirectionalResult:
    """Evaluate one weighted trajectory contribution to Eq. (50)."""

    backward, reweighting = run_backward_trial_path(
        model,
        dividing_surface_state,
        incoming_state=incoming_state,
        temperature=temperature,
        dt=dt,
        approximation=approximation,
        reactant_boundary=reactant_boundary,
        product_boundary=product_boundary,
        max_steps=max_steps,
        rng=backward_rng,
        nuclear_substeps=nuclear_substeps,
        compute_step_energies=compute_step_energies,
        event_convention=event_convention,
        frustrated_hop_policy=frustrated_hop_policy,
        electronic_integrator=electronic_integrator,
    )
    if reweighting is None:
        return BidirectionalResult(
            backward=backward,
            forward=None,
            reweighting=None,
            incoming_state=incoming_state,
            final_state=-1,
            reactive=False,
            forward_crossings=0,
            corrected_weight=0.0,
        )
    forward_initial = TrajectoryState(
        positions=np.asarray(dividing_surface_state.positions).copy(),
        momenta=np.asarray(dividing_surface_state.momenta).copy(),
        active_state=dividing_surface_state.active_state,
        amplitudes=reweighting.dividing_surface_amplitudes.copy(),
    )
    forward = run_forward_rpsh_path(
        model,
        forward_initial,
        temperature=temperature,
        dt=dt,
        approximation=approximation,
        reactant_boundary=reactant_boundary,
        product_boundary=product_boundary,
        max_steps=max_steps,
        rng=forward_rng,
        nuclear_substeps=nuclear_substeps,
        compute_step_energies=compute_step_energies,
        frustrated_hop_policy=frustrated_hop_policy,
        electronic_integrator=electronic_integrator,
    )
    if not forward.completed:
        raise RuntimeError(
            "forward path did not reach either asymptotic boundary; increase max_steps"
        )
    final_centroid = float(forward.centroid_positions[-1])
    reactive = bool(forward.completed and final_centroid >= product_boundary)
    crossings = count_forward_crossings(
        forward.centroid_positions, dividing_surface
    )
    corrected_weight = (
        reweighting.weight / crossings if reactive and crossings > 0 else 0.0
    )
    return BidirectionalResult(
        backward=backward,
        forward=forward,
        reweighting=reweighting,
        incoming_state=incoming_state,
        final_state=forward.final_state.active_state,
        reactive=reactive,
        forward_crossings=crossings,
        corrected_weight=corrected_weight,
    )




