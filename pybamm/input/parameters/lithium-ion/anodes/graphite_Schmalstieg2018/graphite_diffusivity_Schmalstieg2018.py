from pybamm import exp, constants


def graphite_diffusivity_Schmalstieg2018(sto, T):
    """
    Graphite diffusivity as a function of stochiometry [1, 2].

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

    log_D_ref1 = -3.5 * sto - 8.8
    log_D_ref2 = 59.375 * sto ** 3 - 26.563 * sto ** 2 - 8.9125
    log_D_ref3 = -9.7
    # convert cm2s-1 to m2s-1
    D_ref1 = 10 ** log_D_ref1 / 10000
    D_ref2 = 10 ** log_D_ref2 / 10000
    D_ref3 = 10 ** log_D_ref3 / 10000
    
    E_D_s = 2.88e4
    arrhenius = exp(-E_D_s / (constants.R * T)) * exp(E_D_s / (constants.R * 298.15))

    return (sto < 0.2) * (D_ref1 * arrhenius) + (0.2 <= sto < 0.3) * (D_ref2 * arrhenius) + (sto >= 0.3) * (D_ref3 * arrhenius)