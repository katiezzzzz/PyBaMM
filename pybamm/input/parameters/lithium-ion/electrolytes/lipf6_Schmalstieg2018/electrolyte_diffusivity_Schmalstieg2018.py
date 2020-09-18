import pybamm
from pybamm import constants


def electrolyte_diffusivity_Schmalstieg2018(c_e, T):
    """
    Diffusivity of LiPF6 in EC/DMC/EMC as a function of ion concentration [1, 2].

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

    # The diffusivity epends on the electrolyte conductivity
    inputs = {"Electrolyte concentration [mol.m-3]": c_e, "Temperature [K]": T}
    sigma_e = pybamm.FunctionParameter("Electrolyte conductivity [S.m-1]", inputs)

    D_c_e = (sigma_e * constants.R * T) / ((constants.F ** 2) * c_e)

    return D_c_e
