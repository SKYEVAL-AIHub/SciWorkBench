"""Reversible one-dimensional ring-polymer and electronic propagation kernels."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.linalg import expm

from .constants import KB_AU_PER_K as kB

from .electronic import bead_averaged_hamiltonian, centroid_hamiltonian
from .momentum import adjust_ba_momenta, adjust_ca_centroid_momentum
from .rare_event import fssh_probabilities, trial_hop_probabilities


@dataclass(frozen=True)
class TrajectoryState:
    positions: np.ndarray
    momenta: np.ndarray
    active_state: int
    amplitudes: np.ndarray


@dataclass(frozen=True)
class StepResult:
    state: TrajectoryState
    propagator: np.ndarray
    time_derivative_coupling: np.ndarray
    probabilities: np.ndarray
    attempted_state: int | None
    hop_accepted: bool
    energy_before: float
    energy_after: float


def electronic_propagator(
    hamiltonian: np.ndarray, dt: float, *, hermiticity_tolerance: float = 1.0e-10
) -> np.ndarray:
    """Return the unitary constant-H propagator for one midpoint interval."""

    hamiltonian = np.asarray(hamiltonian, dtype=complex)
    if hamiltonian.ndim != 2 or hamiltonian.shape[0] != hamiltonian.shape[1]:
        raise ValueError("hamiltonian must be square")
    if np.linalg.norm(hamiltonian - hamiltonian.conj().T, ord=2) > hermiticity_tolerance:
        raise ValueError("hamiltonian must be Hermitian")
    if hamiltonian.shape == (2, 2):
        center = 0.5 * np.trace(hamiltonian)
        traceless = hamiltonian - center * np.eye(2)
        frequency = float(
            np.sqrt(max(0.0, 0.5 * np.real(np.trace(traceless @ traceless))))
        )
        phase = np.exp(-1j * center * float(dt))
        if frequency < 1.0e-15:
            return phase * np.eye(2, dtype=complex)
        angle = frequency * float(dt)
        return phase * (
            np.cos(angle) * np.eye(2)
            - 1j * np.sin(angle) * traceless / frequency
        )
    return expm(-1j * hamiltonian * float(dt))


def interaction_picture_rk4_propagator(
    energies: np.ndarray,
    time_derivative_couplings: np.ndarray,
    dt: float,
) -> np.ndarray:
    r"""Propagate one interval with the Appendix interaction-picture RK4.

    ``energies`` and ``time_derivative_couplings`` contain the values at the
    beginning, midpoint, and end of the nuclear long-time step.  The diagonal
    phases are integrated from the quadratic interpolant through those three
    values, while the interaction-picture coupling equation is advanced with
    classical fourth-order Runge--Kutta.
    """

    energy = np.asarray(energies, dtype=float)
    coupling = np.asarray(time_derivative_couplings, dtype=complex)
    if energy.ndim != 2 or energy.shape[0] != 3:
        raise ValueError("energies must have shape (3, nstate)")
    nstate = energy.shape[1]
    if coupling.shape != (3, nstate, nstate):
        raise ValueError("time_derivative_couplings must have shape (3, nstate, nstate)")
    if not np.all(np.isfinite(energy)) or not np.all(np.isfinite(coupling)):
        raise ValueError("electronic RK4 inputs must be finite")
    if np.max(np.abs(np.diagonal(coupling, axis1=1, axis2=2))) > 1.0e-12:
        raise ValueError("time-derivative coupling diagonal must vanish")
    if np.max(np.abs(coupling + coupling.swapaxes(1, 2).conj())) > 1.0e-10:
        raise ValueError("time-derivative coupling must be anti-Hermitian")

    step = float(dt)
    theta_mid = step * (5.0 * energy[0] + 8.0 * energy[1] - energy[2]) / 24.0
    theta_end = step * (energy[0] + 4.0 * energy[1] + energy[2]) / 6.0

    def generator(stage_coupling: np.ndarray, theta: np.ndarray) -> np.ndarray:
        phase = np.exp(1j * (theta[:, None] - theta[None, :]))
        return -stage_coupling * phase

    generator_0 = generator(coupling[0], np.zeros(nstate))
    generator_mid = generator(coupling[1], theta_mid)
    generator_end = generator(coupling[2], theta_end)
    identity = np.eye(nstate, dtype=complex)
    k1 = generator_0 @ identity
    k2 = generator_mid @ (identity + 0.5 * step * k1)
    k3 = generator_mid @ (identity + 0.5 * step * k2)
    k4 = generator_end @ (identity + step * k3)
    interaction = identity + step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return np.exp(-1j * theta_end)[:, None] * interaction


@lru_cache(maxsize=None)
def _spring_modes(nbeads: int) -> tuple[np.ndarray, np.ndarray]:
    if nbeads < 1:
        raise ValueError("nbeads must be positive")
    if nbeads == 1:
        # The one-bead (classical) limit has no ring-polymer spring and its
        # only normal mode is the centroid translation.
        return np.zeros(1), np.ones((1, 1))
    laplacian = 2.0 * np.eye(nbeads)
    for bead in range(nbeads):
        laplacian[bead, (bead - 1) % nbeads] -= 1.0
        laplacian[bead, (bead + 1) % nbeads] -= 1.0
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    return np.maximum(eigenvalues, 0.0), eigenvectors


def free_ring_polymer_step(
    positions: np.ndarray,
    momenta: np.ndarray,
    *,
    mass: float,
    temperature: float,
    dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Exactly propagate the free ring polymer, including its internal springs."""

    q = np.asarray(positions, dtype=float)
    p = np.asarray(momenta, dtype=float)
    if q.ndim != 1 or p.shape != q.shape:
        raise ValueError("positions and momenta must be equal 1D arrays")
    if mass <= 0.0 or temperature <= 0.0:
        raise ValueError("mass and temperature must be positive")
    omega_p = q.size * kB * float(temperature)
    indices = np.arange(q.size)
    frequencies = 2.0 * omega_p * np.abs(np.sin(np.pi * indices / q.size))
    q_mode = np.fft.fft(q, norm="ortho")
    p_mode = np.fft.fft(p, norm="ortho")
    phase = frequencies * float(dt)
    cosine = np.cos(phase)
    sine = np.sin(phase)
    q_new = np.empty_like(q_mode, dtype=complex)
    p_new = np.empty_like(p_mode, dtype=complex)
    zero = frequencies < 1.0e-14
    q_new[zero] = q_mode[zero] + float(dt) * p_mode[zero] / float(mass)
    p_new[zero] = p_mode[zero]
    nonzero = ~zero
    q_new[nonzero] = (
        cosine[nonzero] * q_mode[nonzero]
        + sine[nonzero] * p_mode[nonzero] / (mass * frequencies[nonzero])
    )
    p_new[nonzero] = (
        cosine[nonzero] * p_mode[nonzero]
        - mass * frequencies[nonzero] * sine[nonzero] * q_mode[nonzero]
    )
    return (
        np.fft.ifft(q_new, norm="ortho").real,
        np.fft.ifft(p_new, norm="ortho").real,
    )


def ring_polymer_energy(
    model,
    positions: np.ndarray,
    momenta: np.ndarray,
    active_state: int,
    *,
    temperature: float,
) -> float:
    """Return the per-bead active-surface ring-polymer Hamiltonian."""

    q = np.asarray(positions, dtype=float)
    p = np.asarray(momenta, dtype=float)
    if q.ndim != 1 or p.shape != q.shape:
        raise ValueError("positions and momenta must be equal 1D arrays")
    mass = float(model.mass)
    omega_p = q.size * kB * float(temperature)
    spring = 0.5 * mass * omega_p**2 * np.sum((q - np.roll(q, 1)) ** 2)
    potential = np.sum(model.get_energies(q)[:, active_state])
    kinetic = np.sum(p**2) / (2.0 * mass)
    return float((kinetic + spring + potential) / q.size)


def active_surface_step(
    model,
    positions: np.ndarray,
    momenta: np.ndarray,
    active_state: int,
    *,
    temperature: float,
    dt: float,
    substeps: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Symmetric slow-force/free-ring-polymer/slow-force MTS evolution.

    ``substeps`` divides only the fast free-ring-polymer leg, matching the
    multiple-time-step separation stated in the paper's Appendix.  The free
    propagator is analytic here, so changing this value should affect only
    roundoff; keeping the parameter makes the published 10:1 schedule explicit.
    """

    if substeps < 1:
        raise ValueError("substeps must be positive")
    q = np.asarray(positions, dtype=float).copy()
    p = np.asarray(momenta, dtype=float).copy()
    force = model.evaluate(q, need_nac=False, need_wavefun=False).force_diag[
        :, active_state
    ]
    p += 0.5 * float(dt) * force
    small_dt = float(dt) / substeps
    for _ in range(substeps):
        q, p = free_ring_polymer_step(
            q, p, mass=float(model.mass), temperature=temperature, dt=small_dt
        )
    force = model.evaluate(q, need_nac=False, need_wavefun=False).force_diag[
        :, active_state
    ]
    p += 0.5 * float(dt) * force
    return q, p


def effective_hamiltonian(
    model,
    positions: np.ndarray,
    momenta: np.ndarray,
    *,
    approximation: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build the BA or CA Hamiltonian and return its ingredients."""

    q = np.asarray(positions, dtype=float)
    p = np.asarray(momenta, dtype=float)
    velocity = p / float(model.mass)
    approximation = approximation.upper()
    if approximation == "BA":
        snapshot = model.evaluate(q, need_force=False)
        coupling = velocity[:, None, None] * snapshot.nac
        hamiltonian = bead_averaged_hamiltonian(
            snapshot.energies, snapshot.nac, velocity
        )
        return hamiltonian, np.mean(coupling, axis=0), snapshot.energies, snapshot.nac
    if approximation == "CA":
        centroid = float(np.mean(q))
        snapshot = model.evaluate([centroid], need_force=False)
        hamiltonian = centroid_hamiltonian(
            snapshot.energies[0], snapshot.nac[0], float(np.mean(velocity))
        )
        return (
            hamiltonian,
            float(np.mean(velocity)) * snapshot.nac[0],
            model.get_energies(q),
            np.repeat(snapshot.nac, q.size, axis=0),
        )
    raise ValueError("approximation must be 'BA' or 'CA'")


def select_hop(probabilities: np.ndarray, random_number: float) -> int | None:
    """Select one candidate using the usual cumulative-probability rule."""

    probabilities = np.asarray(probabilities, dtype=float)
    if probabilities.ndim != 1 or np.any(probabilities < 0.0):
        raise ValueError("probabilities must be a nonnegative 1D array")
    if probabilities.sum() > 1.0 + 1.0e-12:
        raise ValueError("probability budget exceeds one")
    if not 0.0 <= random_number < 1.0:
        raise ValueError("random_number must lie in [0, 1)")
    cumulative = np.cumsum(probabilities)
    index = int(np.searchsorted(cumulative, random_number, side="right"))
    return None if index >= probabilities.size else index


def _adjust_hop(
    model,
    positions: np.ndarray,
    momenta: np.ndarray,
    old_state: int,
    new_state: int,
    *,
    approximation: str,
) -> tuple[np.ndarray, bool]:
    snapshot = model.evaluate(positions, need_force=False)
    delta_v = float(
        np.mean(snapshot.energies[:, new_state] - snapshot.energies[:, old_state])
    )
    if approximation.upper() == "BA":
        adjustment = adjust_ba_momenta(
            momenta,
            float(model.mass),
            snapshot.nac[:, old_state, new_state],
            delta_v,
        )
    elif approximation.upper() == "CA":
        centroid_nac = model.evaluate(
            [float(np.mean(positions))], need_force=False
        ).nac[0, old_state, new_state]
        adjustment = adjust_ca_centroid_momentum(
            momenta, float(model.mass), centroid_nac, delta_v
        )
    else:
        raise ValueError("approximation must be 'BA' or 'CA'")
    return adjustment.momenta, adjustment.accepted


def _reflect_frustrated_momenta(
    model,
    positions: np.ndarray,
    momenta: np.ndarray,
    old_state: int,
    new_state: int,
    *,
    approximation: str,
) -> np.ndarray:
    """Reverse the NAC-directed momentum component after a frustrated hop."""

    p = np.asarray(momenta, dtype=float).copy()
    approximation = approximation.upper()
    if approximation == "BA":
        direction = model.evaluate(positions, need_force=False).nac[
            :, old_state, new_state
        ]
        norm2 = float(np.dot(direction, direction))
        if norm2 > 1.0e-24:
            p -= 2.0 * float(np.dot(p, direction)) * direction / norm2
        return p
    if approximation == "CA":
        # CA adjusts only the centroid momentum, so reflection reverses the
        # centroid component and leaves all internal normal modes unchanged.
        p -= 2.0 * float(np.mean(p))
        return p
    raise ValueError("approximation must be 'BA' or 'CA'")


def propagate_rpsh_step(
    model,
    state: TrajectoryState,
    *,
    temperature: float,
    dt: float,
    approximation: str,
    random_number: float,
    trial_hop: bool = False,
    frustrated_hop_policy: str = "keep",
    electronic_integrator: str = "midpoint_exponential",
    nuclear_substeps: int = 1,
    compute_energy: bool = True,
) -> StepResult:
    """Propagate one coupled step with either Eq. (53) or ordinary FSSH hops."""

    if frustrated_hop_policy not in {"keep", "reflect"}:
        raise ValueError("frustrated_hop_policy must be 'keep' or 'reflect'")
    if electronic_integrator not in {
        "midpoint_exponential",
        "interaction_picture_rk4",
    }:
        raise ValueError(
            "electronic_integrator must be 'midpoint_exponential' or "
            "'interaction_picture_rk4'"
        )
    q0 = np.asarray(state.positions, dtype=float)
    p0 = np.asarray(state.momenta, dtype=float)
    energy_before = (
        ring_polymer_energy(model, q0, p0, state.active_state, temperature=temperature)
        if compute_energy
        else np.nan
    )
    q1, p1 = active_surface_step(
        model,
        q0,
        p0,
        state.active_state,
        temperature=temperature,
        dt=dt,
        substeps=nuclear_substeps,
    )
    qmid = 0.5 * (q0 + q1)
    pmid = 0.5 * (p0 + p1)
    hamiltonian, coupling, _, _ = effective_hamiltonian(
        model, qmid, pmid, approximation=approximation
    )
    if electronic_integrator == "midpoint_exponential":
        propagator = electronic_propagator(hamiltonian, dt)
    else:
        endpoint_data = [
            effective_hamiltonian(model, q0, p0, approximation=approximation),
            (hamiltonian, coupling, None, None),
            effective_hamiltonian(model, q1, p1, approximation=approximation),
        ]
        stage_hamiltonians = np.asarray([item[0] for item in endpoint_data])
        stage_couplings = np.asarray([item[1] for item in endpoint_data])
        stage_energies = np.real(
            np.diagonal(stage_hamiltonians, axis1=1, axis2=2)
        )
        propagator = interaction_picture_rk4_propagator(
            stage_energies, stage_couplings, dt
        )
    amplitudes = propagator @ np.asarray(state.amplitudes, dtype=complex)
    rho = np.outer(amplitudes, amplitudes.conj())
    interval = abs(float(dt))
    probabilities = (
        trial_hop_probabilities(coupling, state.active_state, interval)
        if trial_hop
        else fssh_probabilities(rho, coupling, state.active_state, interval)
    )
    attempted = select_hop(probabilities, random_number)
    active_state = state.active_state
    accepted = False
    if attempted is not None and attempted != active_state:
        p1, accepted = _adjust_hop(
            model,
            q1,
            p1,
            active_state,
            attempted,
            approximation=approximation,
        )
        if accepted:
            active_state = attempted
        elif frustrated_hop_policy == "reflect":
            p1 = _reflect_frustrated_momenta(
                model,
                q1,
                p1,
                active_state,
                attempted,
                approximation=approximation,
            )
    energy_after = (
        ring_polymer_energy(model, q1, p1, active_state, temperature=temperature)
        if compute_energy
        else np.nan
    )
    return StepResult(
        state=TrajectoryState(q1, p1, active_state, amplitudes),
        propagator=propagator,
        time_derivative_coupling=coupling,
        probabilities=probabilities,
        attempted_state=attempted,
        hop_accepted=accepted,
        energy_before=energy_before,
        energy_after=energy_after,
    )




