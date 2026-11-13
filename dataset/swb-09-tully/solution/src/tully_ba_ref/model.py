"""Self-contained implementation of Eq. (57) from Shushkov, Li, and Tully.

All quantities are in atomic units.  The globally smooth analytic adiabatic
gauge avoids numerical eigenvector sign flips and therefore makes the bead
average in Eq. (36) deterministic.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ModelSnapshot:
    x: np.ndarray
    energies: np.ndarray | None = None
    wavefun: np.ndarray | None = None
    force: np.ndarray | None = None
    force_diag: np.ndarray | None = None
    nac: np.ndarray | None = None


class Tully2012ReactionModel:
    """One-dimensional, two-state reaction model in the adiabatic basis."""

    def __init__(self, representation: str = "adiabatic"):
        if representation != "adiabatic":
            raise ValueError("the standalone reference supports the adiabatic basis only")
        self.representation = representation
        self.nstate = 2
        self.mass = 2000.0
        self.A = 5.0e-3
        self.B = 2.0
        self.C = 1.0e-5
        self.D = -2.8
        self.E = 1.2
        self.F = 1.0e-3
        self.G = 2.0

    @staticmethod
    def _points(x) -> tuple[np.ndarray, bool]:
        array = np.asarray(x, dtype=float)
        scalar = array.ndim == 0
        return np.atleast_1d(array), scalar

    def V(self, x) -> np.ndarray:
        points, scalar = self._points(x)
        v11 = self.A * np.exp(-self.B * points**2)
        v22 = self.C * np.exp(2.0 * self.E * (points - self.D))
        v22 += 2.0 * self.C * np.exp(self.E * (points - self.D))
        v12 = self.F * np.exp(-self.G * points**2)
        out = np.empty((points.size, 2, 2), dtype=float)
        out[:, 0, 0] = v11
        out[:, 1, 1] = v22
        out[:, 0, 1] = out[:, 1, 0] = v12
        return out[0] if scalar else out

    def dV(self, x) -> np.ndarray:
        points, scalar = self._points(x)
        v11 = -2.0 * self.A * self.B * points * np.exp(-self.B * points**2)
        v22 = 2.0 * self.E * self.C * np.exp(2.0 * self.E * (points - self.D))
        v22 += 2.0 * self.C * self.E * np.exp(self.E * (points - self.D))
        v12 = -2.0 * self.F * self.G * points * np.exp(-self.G * points**2)
        out = np.empty((points.size, 2, 2), dtype=float)
        out[:, 0, 0] = v11
        out[:, 1, 1] = v22
        out[:, 0, 1] = out[:, 1, 0] = v12
        return out[0] if scalar else out

    def evaluate(
        self,
        x,
        *,
        wavefunc_for_nac=None,
        energies_for_nac=None,
        need_force: bool = True,
        need_nac: bool = True,
        need_wavefun: bool = True,
    ) -> ModelSnapshot:
        del wavefunc_for_nac, energies_for_nac
        points, _ = self._points(x)
        diabatic = self.V(points)
        derivative = self.dV(points)
        v11 = diabatic[:, 0, 0]
        v22 = diabatic[:, 1, 1]
        v12 = diabatic[:, 0, 1]
        delta = v11 - v22
        radius = np.sqrt(0.25 * delta**2 + v12**2)
        mean = 0.5 * (v11 + v22)
        energies = np.column_stack((mean - radius, mean + radius))

        theta = 0.5 * np.arctan2(2.0 * v12, delta)
        sine = np.sin(theta)
        cosine = np.cos(theta)
        wavefun = np.empty((points.size, 2, 2), dtype=float)
        wavefun[:, :, 0] = np.column_stack((-sine, cosine))
        wavefun[:, :, 1] = np.column_stack((cosine, sine))

        force = None
        force_diag = None
        if need_force:
            force = -np.einsum(
                "bpi,bpq,bqj->bij", wavefun, derivative, wavefun, optimize=True
            )
            force_diag = np.diagonal(force, axis1=1, axis2=2)

        nac = None
        if need_nac:
            delta_prime = derivative[:, 0, 0] - derivative[:, 1, 1]
            coupling_prime = derivative[:, 0, 1]
            denominator = delta**2 + 4.0 * v12**2
            theta_prime = (delta * coupling_prime - v12 * delta_prime) / denominator
            nac = np.zeros((points.size, 2, 2), dtype=float)
            nac[:, 0, 1] = theta_prime
            nac[:, 1, 0] = -theta_prime

        return ModelSnapshot(
            x=points,
            energies=energies,
            wavefun=wavefun if need_wavefun else None,
            force=force,
            force_diag=force_diag,
            nac=nac,
        )

    def get_energies(self, x) -> np.ndarray:
        return self.evaluate(x, need_force=False, need_nac=False, need_wavefun=False).energies

    def get_force_diag(self, x) -> np.ndarray:
        return self.evaluate(x, need_force=True, need_nac=False, need_wavefun=False).force_diag


class Tully2012AdiabaticNoHopModel(Tully2012ReactionModel):
    """Same adiabatic surfaces with derivative couplings disabled."""

    def evaluate(self, x, **kwargs) -> ModelSnapshot:
        snapshot = super().evaluate(x, **kwargs)
        if snapshot.nac is None:
            return snapshot
        return ModelSnapshot(
            x=snapshot.x,
            energies=snapshot.energies,
            wavefun=snapshot.wavefun,
            force=snapshot.force,
            force_diag=snapshot.force_diag,
            nac=np.zeros_like(snapshot.nac),
        )



