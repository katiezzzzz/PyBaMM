from pybamm import exp, constants, Parameter


def graphite_negative_particle_distribution_in_x_Ecker2015(x):
    """
    Negative particle distribution in x

    References
    ----------
    .. [1] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery i. determination of parameters." Journal of the
    Electrochemical Society 162.9 (2015): A1836-A1848.
    .. [2] Ecker, Madeleine, et al. "Parameterization of a physico-chemical model of
    a lithium-ion battery ii. model validation." Journal of The Electrochemical
    Society 162.9 (2015): A1849-A1857.
    .. [3] Richardson, Giles, et. al. "Generalised single particle models for
    high-rate operation of graded lithium-ion electrodes: Systematic derivation
    and validation." Electrochemica Acta 339 (2020): 135862

    Parameters
    ----------
    x : :class:`pybamm.Symbol`
        Through-cell distance (x_n) [m]

    Returns
    -------
    :class:`pybamm.Symbol`
        Negative particle distribution in x
    """
    
    NegativeElectrodeThickness = 7.4e-05
    
    return (x / NegativeElectrodeThickness) + 1/2