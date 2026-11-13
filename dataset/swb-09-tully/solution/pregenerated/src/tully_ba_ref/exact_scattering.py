"""Stationary coupled-channel quantum scattering on a finite-difference grid."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from .constants import KB_AU_PER_K as kB


@dataclass(frozen=True)
class ThermalScatteringRate:
    temperature: float
    energies: np.ndarray
    quadrature_weights: np.ndarray
    transmission_by_incoming_state: np.ndarray
    rate_by_incoming_state: np.ndarray
    total_rate: float
    grid: np.ndarray


def _retarded_lead_self_energy(energy: float, potential: float, hopping: float) -> complex:
    """Surface self-energy for a semi-infinite nearest-neighbour 1D lead."""

    z = complex(energy - (2.0 * hopping + potential), 1.0e-14)
    root = np.lib.scimath.sqrt(z * z - 4.0 * hopping * hopping)
    candidates = (0.5 * (z - root), 0.5 * (z + root))
    physical = min(candidates, key=lambda value: (abs(value), value.imag > 0.0))
    if physical.imag > 0.0:
        physical = physical.conjugate()
    return physical


def _lead_matrices(
    energy: float, potential_matrix: np.ndarray, hopping: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    eigenvalues, eigenvectors = np.linalg.eigh(potential_matrix)
    sigma_channels = np.array(
        [_retarded_lead_self_energy(energy, value, hopping) for value in eigenvalues]
    )
    gamma_channels = np.maximum(0.0, -2.0 * np.imag(sigma_channels))
    sigma = (eigenvectors * sigma_channels) @ eigenvectors.conj().T
    gamma = (eigenvectors * gamma_channels) @ eigenvectors.conj().T
    return sigma, gamma, gamma_channels, eigenvectors


def state_resolved_transmission_fd(
    potential,
    *,
    energy: float,
    mass: float,
    xmin: float,
    xmax: float,
    ngrid: int,
) -> np.ndarray:
    """Compute left-incoming transmission probabilities with open lead self-energies."""

    if energy < 0.0 or mass <= 0.0 or ngrid < 5 or xmin >= xmax:
        raise ValueError("invalid scattering-grid parameters")
    grid = np.linspace(xmin, xmax, ngrid)
    dx = float(grid[1] - grid[0])
    hopping = 1.0 / (2.0 * mass * dx * dx)
    matrices = np.asarray(potential(grid), dtype=complex)
    if matrices.ndim != 3 or matrices.shape[0] != ngrid or matrices.shape[1] != matrices.shape[2]:
        raise ValueError("potential must return (ngrid, nstate, nstate)")
    nstate = matrices.shape[1]
    kinetic = sparse.diags(
        [-hopping * np.ones(ngrid - 1), 2.0 * hopping * np.ones(ngrid), -hopping * np.ones(ngrid - 1)],
        offsets=(-1, 0, 1),
        format="csc",
    )
    hamiltonian = sparse.kron(kinetic, sparse.eye(nstate), format="csc")
    hamiltonian += sparse.block_diag(list(matrices), format="csc")
    sigma_left, gamma_left, gamma_left_channels, left_vectors = _lead_matrices(
        energy, matrices[0], hopping
    )
    sigma_right, gamma_right, _, _ = _lead_matrices(
        energy, matrices[-1], hopping
    )
    size = ngrid * nstate
    operator = sparse.eye(size, dtype=complex, format="lil") * complex(energy, 1.0e-12)
    operator -= hamiltonian
    operator[:nstate, :nstate] -= sigma_left
    operator[-nstate:, -nstate:] -= sigma_right
    factor = splu(operator.tocsc())
    sources = np.zeros((size, nstate), dtype=complex)
    sources[:nstate, :] = left_vectors
    green_right_left = factor.solve(sources)[-nstate:, :]
    transmission = np.empty(nstate)
    for incoming in range(nstate):
        outgoing = green_right_left[:, incoming]
        transmission[incoming] = gamma_left_channels[incoming] * float(
            np.real(np.vdot(outgoing, gamma_right @ outgoing))
        )
    transmission[np.abs(transmission) < 1.0e-12] = 0.0
    return np.maximum(transmission, 0.0)


def thermal_nonadiabatic_rate(
    model,
    *,
    temperature: float,
    xmin: float,
    xmax: float,
    ngrid: int,
    nquadrature: int,
) -> ThermalScatteringRate:
    r"""Thermally integrate exact state-resolved transmission probabilities.

    With \(x=\beta E\), the dimensionless transmission average is evaluated
    directly by Gauss-Laguerre quadrature,
    \(\int_0^\infty e^{-x}T_i(x/\beta)\,dx\).
    """

    if temperature <= 0.0 or nquadrature < 4:
        raise ValueError("temperature and quadrature order must be positive")
    beta = 1.0 / (kB * float(temperature))
    nodes, weights = np.polynomial.laguerre.laggauss(nquadrature)
    energies = nodes / beta
    transmissions = np.array(
        [
            state_resolved_transmission_fd(
                model.V,
                energy=float(energy),
                mass=float(model.mass),
                xmin=xmin,
                xmax=xmax,
                ngrid=ngrid,
            )
            for energy in energies
        ]
    )
    free_flux = 0.5 / np.sqrt(2.0 * np.pi * float(model.mass) * beta)
    rates = free_flux * np.sum(weights[:, None] * transmissions, axis=0)
    return ThermalScatteringRate(
        temperature=float(temperature),
        energies=energies,
        quadrature_weights=weights,
        transmission_by_incoming_state=transmissions,
        rate_by_incoming_state=rates,
        total_rate=float(np.sum(rates)),
        grid=np.linspace(xmin, xmax, ngrid),
    )


def thermal_adiabatic_rate(
    model,
    *,
    temperature: float,
    xmin: float,
    xmax: float,
    ngrid: int,
    nquadrature: int,
) -> ThermalScatteringRate:
    """Thermal quantum rate on the lower adiabatic PES without derivative couplings."""

    beta = 1.0 / (kB * float(temperature))
    nodes, weights = np.polynomial.laguerre.laggauss(nquadrature)
    energies = nodes / beta

    def lower_surface(grid):
        lower = np.linalg.eigvalsh(model.V(grid))[:, 0]
        return lower[:, None, None]

    transmissions = np.array(
        [
            state_resolved_transmission_fd(
                lower_surface,
                energy=float(energy),
                mass=float(model.mass),
                xmin=xmin,
                xmax=xmax,
                ngrid=ngrid,
            )
            for energy in energies
        ]
    )
    free_flux = 0.5 / np.sqrt(2.0 * np.pi * float(model.mass) * beta)
    rates = free_flux * np.sum(weights[:, None] * transmissions, axis=0)
    return ThermalScatteringRate(
        temperature=float(temperature),
        energies=energies,
        quadrature_weights=weights,
        transmission_by_incoming_state=transmissions,
        rate_by_incoming_state=rates,
        total_rate=float(np.sum(rates)),
        grid=np.linspace(xmin, xmax, ngrid),
    )




