from pybamm import exp, constants


def nmc_diffusivity_Schmalstieg2018(sto, T):
    """
    NMC diffusivity as a function of stoichiometry [1, 2].

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
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Solid diffusivity
    """

    log_D_ref = -1.682 * sto - 9.127
    # convert cm2s-1 to m2s-1
    D_ref = 10 ** log_D_ref / 10000
    E_D_s = 4.89e4
    arrhenius = exp(-E_D_s / (constants.R * T)) * exp(E_D_s / (constants.R * 298.15))

    return D_ref * arrhenius
