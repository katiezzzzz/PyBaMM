#
# Tests for the lithium-ion DFN model
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm
import tests

import numpy as np
import unittest


class TestDFN(unittest.TestCase):
    def test_basic_processing(self):
        model = pybamm.lithium_ion.DFN()
        var = pybamm.standard_spatial_vars
        model.default_var_pts = {
            var.x_n: 3,
            var.x_s: 3,
            var.x_p: 3,
            var.r_n: 3,
            var.r_p: 3,
        }

        modeltest = tests.StandardModelTest(model)
        modeltest.test_all()

    def test_optimisations(self):
        model = pybamm.lithium_ion.DFN()
        optimtest = tests.OptimisationsTest(model)

        original = optimtest.evaluate_model()
        simplified = optimtest.evaluate_model(simplify=True)
        using_known_evals = optimtest.evaluate_model(use_known_evals=True)
        simp_and_known = optimtest.evaluate_model(simplify=True, use_known_evals=True)
        np.testing.assert_array_almost_equal(original, simplified)
        np.testing.assert_array_almost_equal(original, using_known_evals)
        np.testing.assert_array_almost_equal(original, simp_and_known)


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()