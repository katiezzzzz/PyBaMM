from pybamm import exp, constants, Parameter


def nmc_electrolyte_exchange_current_density_Schmalstieg2018(c_e, c_s_surf, T):
    """
    Exchange-current density for Butler-Volmer reactions between NMC and LiPF6 in
    EC/DMC/EMC.

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
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    
    k_ref = 5.03 / (constants.F * 1000 ** 0.5 * 24195 ** 0.5 * (48390 - 24195) ** 0.5)
    m_ref = constants.F * k_ref  # (A/m2)(mol/m3)**1.5 - includes ref concentrations
    
    E_r = 78100
    arrhenius = exp(-E_r / (constants.R * T)) * exp(E_r / (constants.R * 298.15))

    c_p_max = Parameter("Maximum concentration in positive electrode [mol.m-3]")

    return (
        m_ref * arrhenius * c_e ** 0.5 * c_s_surf ** 0.5 * (c_p_max - c_s_surf) ** 0.5
    )
