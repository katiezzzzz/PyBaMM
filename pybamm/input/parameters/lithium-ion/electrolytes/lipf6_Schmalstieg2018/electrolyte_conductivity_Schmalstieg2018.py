from pybamm import exp, constants


def electrolyte_conductivity_Schmalstieg2018(c_e, T):
    """
    Conductivity of LiPF6 in EC/DMC/EMC as a function of ion concentration [1, 2].

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
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    # mol/m^3 to mol/l
    cm = 1e-3 * c_e

    # regression function
    regression = -5.384 + 0.03213 * T - 0.00368 * cm * T + 1.320 * exp(-2.235 * cm)

    # calculate conducvitity
    sigma_e = (regression ** 2) * cm

    return sigma_e
