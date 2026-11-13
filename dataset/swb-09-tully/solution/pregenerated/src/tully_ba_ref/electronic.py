"""Electronic Hamiltonians defining the centroid and bead approximations."""

from __future__ import annotations

import numpy as np


def _validate_bead_inputs(
    energies: np.ndarray, nac: np.ndarray, velocities: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    energies = np.asarray(energies, dtype=float)
    nac = np.asarray(nac, dtype=float)
    velocities = np.asarray(velocities, dtype=float)
    if energies.ndim != 2:
        raise ValueError("energies must have shape (nbead, nstate)")
    nbead, nstate = energies.shape
    if nac.shape != (nbead, nstate, nstate):
        raise ValueError("nac must have shape (nbead, nstate, nstate)")
    if velocities.shape != (nbead,):
        raise ValueError("velocities must have shape (nbead,)")
    return energies, nac, velocities


def bead_averaged_hamiltonian(
    energies: np.ndarray, nac: np.ndarray, velocities: np.ndarray
) -> np.ndarray:
    r"""Return the BA electronic Hamiltonian of Eq. (36).

    In one dimension,
    \(H^{\mathrm{BA}}_{ij}=\delta_{ij}\langle E_i\rangle_P
    -i\langle \dot R d_{ij}\rangle_P\).
    The caller must supply bead NACs in one common electronic gauge.
    """

    energies, nac, velocities = _validate_bead_inputs(energies, nac, velocities)
    hamiltonian = np.diag(np.mean(energies, axis=0)).astype(complex)
    hamiltonian -= 1j * np.mean(velocities[:, None, None] * nac, axis=0)
    return hamiltonian


def centroid_hamiltonian(
    energies: np.ndarray, nac: np.ndarray, centroid_velocity: float
) -> np.ndarray:
    r"""Return the CA electronic Hamiltonian of Eq. (35)."""

    energies = np.asarray(energies, dtype=float)
    nac = np.asarray(nac, dtype=float)
    if energies.ndim != 1:
        raise ValueError("centroid energies must have shape (nstate,)")
    if nac.shape != (energies.size, energies.size):
        raise ValueError("centroid nac must have shape (nstate, nstate)")
    hamiltonian = np.diag(energies).astype(complex)
    hamiltonian -= 1j * float(centroid_velocity) * nac
    return hamiltonian


def mean_active_surface_gap(
    bead_energies: np.ndarray, old_state: int, new_state: int
) -> float:
    """Mean ring-polymer potential change for a whole-ring hop."""

    bead_energies = np.asarray(bead_energies, dtype=float)
    if bead_energies.ndim != 2:
        raise ValueError("bead_energies must have shape (nbead, nstate)")
    return float(np.mean(bead_energies[:, new_state] - bead_energies[:, old_state]))




