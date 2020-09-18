from pybamm import tanh, exp


def nmc_ocp_Schmalstieg2018_function(sto):
    """
    NMC OCP as a function of stochiometry [1, 2].

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
    sto : :class:`pybamm.Symbol`
       Stochiometry of material (li-fraction)

    """
    
    a = 2.221
    b = 1.696
    c = 0.1578
    d = 153.6
    e = 3.013

    u_eq = (a * exp(-b * sto) + c) * tanh(d * (1 - sto)) + e
    
    return u_eq
