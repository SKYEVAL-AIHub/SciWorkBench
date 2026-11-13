from __future__ import annotations

import unittest

import numpy as np

from tully_ba_ref.electronic import bead_averaged_hamiltonian, centroid_hamiltonian
from tully_ba_ref.model import Tully2012ReactionModel
from tully_ba_ref.momentum import adjust_ba_momenta, mean_kinetic_energy
from tully_ba_ref.rare_event import backward_step_weight


class StandaloneReferenceTests(unittest.TestCase):
    def setUp(self):
        self.model = Tully2012ReactionModel()

    def test_eq57_potential(self):
        value = self.model.V(0.0)
        self.assertAlmostEqual(value[0, 0], 0.005)
        self.assertAlmostEqual(value[0, 1], 0.001)
        expected = 1e-5 * np.exp(2 * 1.2 * 2.8) + 2e-5 * np.exp(1.2 * 2.8)
        self.assertAlmostEqual(value[1, 1], expected)

    def test_analytic_gauge(self):
        points = np.linspace(-3.0, 2.0, 101)
        snapshot = self.model.evaluate(points)
        np.testing.assert_allclose(snapshot.energies, np.linalg.eigvalsh(self.model.V(points)))
        np.testing.assert_allclose(snapshot.nac + snapshot.nac.transpose(0, 2, 1), 0.0, atol=1e-13)

    def test_one_bead_ba_equals_ca(self):
        snapshot = self.model.evaluate([0.1])
        velocity = np.array([0.003])
        ba = bead_averaged_hamiltonian(snapshot.energies, snapshot.nac, velocity)
        ca = centroid_hamiltonian(snapshot.energies[0], snapshot.nac[0], velocity[0])
        np.testing.assert_allclose(ba, ca)

    def test_ba_hop_conserves_mean_energy(self):
        momenta = np.array([3.0, 2.0, 4.0, 2.5])
        direction = np.array([1.0, -0.5, 0.75, 1.25])
        delta_v = -1e-3
        before = mean_kinetic_energy(momenta, 2000.0)
        result = adjust_ba_momenta(momenta, 2000.0, direction, delta_v)
        self.assertTrue(result.accepted)
        after = mean_kinetic_energy(result.momenta, 2000.0)
        self.assertAlmostEqual(after - before + delta_v, 0.0, places=12)

    def test_change_of_measure_identity(self):
        g = np.array([0.0, 0.2])
        f = np.array([0.0, 0.1])
        expectation = f[1] * backward_step_weight(g, f, 1)
        expectation += (1.0 - f.sum()) * backward_step_weight(g, f, None)
        self.assertAlmostEqual(expectation, 1.0)


if __name__ == "__main__":
    unittest.main()



