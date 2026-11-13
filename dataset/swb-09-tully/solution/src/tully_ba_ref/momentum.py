"""Energy-conserving momentum changes for whole-ring BA and CA hops."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MomentumAdjustment:
    momenta: np.ndarray
    accepted: bool
    coefficient: float | None
    discriminant: float
    energy_error: float
    reason: str


def mean_kinetic_energy(momenta: np.ndarray, mass: float) -> float:
    """Ring-polymer kinetic energy in the paper's per-bead convention."""

    momenta = np.asarray(momenta, dtype=float)
    return float(np.mean(momenta**2) / (2.0 * float(mass)))


def _smallest_root(a: float, b: float, c: float, tol: float) -> tuple[float | None, float]:
    disc = b * b - 4.0 * a * c
    scale = max(b * b, abs(4.0 * a * c), 1.0)
    if disc < -tol * scale:
        return None, disc
    disc = max(disc, 0.0)
    sqrt_disc = float(np.sqrt(disc))
    roots = ((-b + sqrt_disc) / (2.0 * a), (-b - sqrt_disc) / (2.0 * a))
    return min(roots, key=abs), disc


def adjust_ba_momenta(
    momenta: np.ndarray,
    mass: float,
    bead_nac_direction: np.ndarray,
    delta_v_mean: float,
    *,
    tol: float = 1.0e-12,
) -> MomentumAdjustment:
    r"""Adjust all bead momenta parallel to the multi-bead NAC vector.

    The update is \(p'_k=p_k+\gamma d_k\).  It conserves
    \(P^{-1}\sum_k p_k^2/(2M)+\bar V_s\), consistent with a
    whole-ring hop on the bead-averaged active surface.
    """

    p = np.asarray(momenta, dtype=float)
    direction = np.asarray(bead_nac_direction, dtype=float)
    if p.ndim != 1 or direction.shape != p.shape:
        raise ValueError("momenta and bead_nac_direction must be equal 1D arrays")
    mass = float(mass)
    if mass <= 0.0:
        raise ValueError("mass must be positive")
    c = float(delta_v_mean)
    norm2 = float(np.dot(direction, direction))
    if norm2 <= tol:
        if abs(c) <= tol:
            return MomentumAdjustment(p.copy(), True, 0.0, 0.0, 0.0, "zero gap")
        return MomentumAdjustment(p.copy(), False, None, -np.inf, c, "zero NAC direction")

    nbead = p.size
    a = norm2 / (2.0 * mass * nbead)
    b = float(np.dot(p, direction)) / (mass * nbead)
    gamma, disc = _smallest_root(a, b, c, tol)
    if gamma is None:
        return MomentumAdjustment(p.copy(), False, None, disc, c, "frustrated hop")
    adjusted = p + gamma * direction
    error = mean_kinetic_energy(adjusted, mass) - mean_kinetic_energy(p, mass) + c
    return MomentumAdjustment(adjusted, True, gamma, disc, error, "accepted")


def adjust_ca_centroid_momentum(
    momenta: np.ndarray,
    mass: float,
    centroid_nac_direction: float,
    delta_v_mean: float,
    *,
    tol: float = 1.0e-12,
) -> MomentumAdjustment:
    r"""Adjust only the centroid normal-mode momentum for a CA hop.

    In this one-dimensional implementation every bead receives the same
    increment, so all internal-mode momenta remain exactly unchanged.
    """

    p = np.asarray(momenta, dtype=float)
    if p.ndim != 1:
        raise ValueError("momenta must be a 1D array")
    mass = float(mass)
    if mass <= 0.0:
        raise ValueError("mass must be positive")
    direction = float(centroid_nac_direction)
    c = float(delta_v_mean)
    if abs(direction) <= tol:
        if abs(c) <= tol:
            return MomentumAdjustment(p.copy(), True, 0.0, 0.0, 0.0, "zero gap")
        return MomentumAdjustment(p.copy(), False, None, -np.inf, c, "zero NAC direction")

    pbar = float(np.mean(p))
    a = direction * direction / (2.0 * mass)
    b = pbar * direction / mass
    gamma, disc = _smallest_root(a, b, c, tol)
    if gamma is None:
        return MomentumAdjustment(p.copy(), False, None, disc, c, "frustrated hop")
    adjusted = p + gamma * direction
    error = mean_kinetic_energy(adjusted, mass) - mean_kinetic_energy(p, mass) + c
    return MomentumAdjustment(adjusted, True, gamma, disc, error, "accepted")




