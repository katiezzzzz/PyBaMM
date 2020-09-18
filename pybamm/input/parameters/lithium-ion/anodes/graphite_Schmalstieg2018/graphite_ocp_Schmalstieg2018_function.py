from pybamm import exp, tanh


def graphite_ocp_Schmalstieg2018_function(sto):
    """
    Graphite OCP as a function of stochiometry [1, 2].

    References
    ----------
    .. [1] Schmalstieg, Johannes, et al. "Full cell parameterization of a high-power lithium-ion 
    battery for a physico-chemical model: part i. physical and electrochemical parameters." 
    Journal of the Electrochemical Society 165.16 (2018): A3799-A3810.
    .. [2] Schmalstieg, Johannes, et al. "Full cell parameterization of a high-power lithium-ion 
    battery for a physico-chemical model: part ii. thermal parameters and validation." 
    Journal of The Electrochemical Society 165.16 (2018): A3811-A3819.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open circuit potential
    """
    A1 = 1.749
    A2 = 0.607
    A3 = 0.01902
    A4 = 0.07714
    k1 = 99.05
    k2 = 15.72
    k3 = 33.91
    k4 = 29.74
    B = 0.24
    C = 0.5237
    D = 0.03215
    
    correction = A2 * exp(-k2 * 0.125) + A3 * tanh(k3 * (C - 0.125)) + A4 * tanh(k4 * (1 - 0.125)) + D - (A1 * exp(-k1 * 0.125) + B)
    u_eq1 = A1 * exp(-k1 * sto) + B + correction
    u_eq2 = A2 * exp(-k2 * sto) + A3 * tanh(k3 * (C - sto)) + A4 * tanh(k4 * (1 - sto)) + D

    return (sto < 0.125) * u_eq1 + (sto >= 0.125) * u_eq2
