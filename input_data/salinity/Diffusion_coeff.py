import numpy as np


def calculate_viscosity_np(C, T):
    """
    taken from Batzle_Wang paper
    :param C:salinity in ppm
    :param T:Temperature in degree celsius
    :return:viscosity in cP
    """

    viscosity = 0.1 + (0.333 * C) + (1.65 + 91.9 * C ** 3) * np.exp(-(0.42 * (C ** 0.8 - 0.17) ** 2 + 0.045) * T ** 0.8)

    return viscosity


def calculate_diffusion_coeff(T_cal):
    """
     taken from paper written by (Simpson and Carr, 1958)
    :param T_cal:Temperature in degree
    :param D_ref :diffusion coefficient in m^2/s at 25 degrees celsius
    :param viscosity_ref: viscosity of water at 25 degrees celsius in cP
    :param T_ref:temperature of reference in kelvin
    :param viscosity_cal:calculated viscosity in cp using Batzle_wang
    :return:D_cal:calculated diffusion coefficient in the new temperature
    """
    T_cal += 273.15
    D_ref = 2.03 * 10 ** -9
    T_ref = 298.15
    viscosity_ref = 0.890
    viscosity_cal = calculate_viscosity_np(0, T_cal)
    D_cal = ((D_ref * viscosity_ref) / T_ref) * (T_cal / viscosity_cal)

    return D_cal


