from pybamm import exp, constants


def graphite_diffusivity_Ecker2015_coin(sto, T):
    """
    Graphite diffusivity as a function of stochiometry [1, 2].

    References
    ----------
     .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.

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
   # "piecewise coin" function for coin model
    y0 = -15.03
    A4 = 1.187
    sigma4 = 0.02648

    log_D_ref = y0 + (
        A4 * exp(-((sto - 1) / sigma4) ** 2)
    )
    
    limit = 0.092
    A1 = 14.4e-13
    A2 = 17.5e-13
    B1 = -2.3
    B2 = -16.94
    C = 0.83e-15
    
    D_ref_low = A1 * exp(B1 * sto) + A2 * exp(B2 * limit) + C - A1 * exp(B1 * limit) + 10 ** log_D_ref
    D_ref_high = A2 * exp(B2 * sto) + C + 10 ** log_D_ref
    E_D_s = 3.03e4
    arrhenius = exp(-E_D_s / (constants.R * T)) * exp(E_D_s / (constants.R * T))
    
    ratio = 2.1e-05 / 1.37e-05

    return (sto < limit) * (D_ref_low * arrhenius * ((1.5 * ratio) ** 2)) + (sto >= limit) * (D_ref_high * arrhenius * ((1.5 * ratio) ** 2))