#
# Compare basic models with full models
#
import pybamm

import numpy as np
import unittest


class TestCompareBasicModels(unittest.TestCase):
    def test_compare_full(self):
        basic_full = pybamm.lead_acid.BasicFull()
        full = pybamm.lead_acid.Full({"convection": "uniform transverse"})

        parameter_values = pybamm.ParameterValues(
            chemistry=pybamm.parameter_sets.Sulzer2019
        )
        parameter_values["Current function [A]"] = 10

        # Solve basic Full mode
        basic_sim = pybamm.Simulation(
            basic_full, solver=pybamm.CasadiSolver(), parameter_values=parameter_values,
        )
        t_eval = np.linspace(0, 400)
        basic_sim.solve(t_eval)
        basic_sol = basic_sim.solution

        # Solve main Full model
        sim = pybamm.Simulation(
            full, solver=pybamm.CasadiSolver(), parameter_values=parameter_values
        )
        t_eval = np.linspace(0, 400)
        sim.solve(t_eval)
        sol = sim.solution

        # Compare solution data
        # np.testing.assert_array_almost_equal(basic_sol.y, sol.y, decimal=4)
        np.testing.assert_array_almost_equal(basic_sol.t, sol.t, decimal=4)
        # Compare variables
        for name in basic_full.variables:
            np.testing.assert_array_almost_equal(
                basic_sol[name].entries, sol[name].entries, decimal=4
            )


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()