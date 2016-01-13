__author__ = 'elco'

import numpy as np
import pdb


def He_diffusion_Meesters_and_Dunai_2002(t, D, radius, Ur0, U_function='constant', shape='sphere',
                                         decay_constant=1e-8, n_eigenmodes=15, x=0, all_timesteps=True,
                                         alpha_ejection=True, stopping_distance=20e-6):

    """
    Helium diffusion model by Meesters and Dunai (2002)

    t : numpy array
        time (s)
    D : numpy array
        diffusivity for each timestep (m2 s-1)
    radius : float
        radius of diffusion domain
    U_function : string, optional
        'constant' or 'exponential'
    shape : string, optional
        shape of the modeled diffusion domain, default is 'sphere'
    n_eigenmodes : int
        number of eigenmodes to evaluate
    x : float
        initial concentration of daughter isotope

    """

    # find shape params
    n = np.arange(1, n_eigenmodes+1, 1)
    if shape == 'sphere':
        mu = (n * np.pi / radius) ** 2
        gamma = 6.0 / (n * np.pi) ** 2
    elif shape == 'finite cylinder':
        print 'need to figure out bessel functions in numpy'
        # http://docs.scipy.org/doc/scipy-0.13.0/reference/special.html
    elif shape == 'infinite cylinder':
        print 'need to figure out bessel functions in numpy'

    # experimental, implement alpha ejection algorithm
    # this is bascially adjusting the gamma term, see eq. 24 in Meesters and Dunai (2002) pt. 2
    # sigma is alpha ejection distance
    # acc. Farley et al. (1996) Geochem Cosm. Acta, varies between 10-30 um for apatites
    if alpha_ejection is True:
        a = radius
        k = n * np.pi / a
        sigma = stopping_distance
        gamma_old = gamma
        gamma = 3.0 / (n * np.pi)**2 * (1.0 - sigma / (2 * a) + 1 / (n * np.pi) * 1 / (k * sigma) * (1.0 - np.cos(k * sigma)) + 1.0 / (k * sigma) * np.sin(k * sigma))

        pdb.set_trace()

    nt = len(t)

    # calculate decay time
    decay_time = 1.0 / decay_constant

    # calculate F function (eq ..)
    if U_function == 'constant':
        F = t
    elif U_function == 'exponential':
        F = decay_time * (1.0 - np.exp(-t / decay_time))
    else:
        raise ValueError('please supply value for U_function, choose either "constant" or "exponential"')

    # eq. 5, time integration of diffusivity:
    xi = np.zeros(nt)
    for j in xrange(1, nt):
        xi[j] = xi[j-1] + (D[j-1] + D[j]) / 2.0 * (t[j] - t[j-1])

    # eq. 6
    Fa = (F[1:] - F[:-1]) / (xi[1:] - xi[:-1])

    # iterate over all eigenmodes:
    beta = np.zeros((n_eigenmodes, nt))

    if all_timesteps is True:

        cn = np.zeros((nt, n_eigenmodes))

        for J in xrange(2, nt):

            beta[:, :] = 0
            #beta_sum[:] = 0

            for n in xrange(n_eigenmodes):

                # eq. 8
                beta[n, :J] = np.exp(-mu[n] * (xi[J] - xi[:J]))

                # right hand side of eq. 7:
                beta_sum = np.sum((beta[n, 1:J] - beta[n, :(J-1)]) * Fa[:(J-1)])

                # eq. 7
                cn[J, n] = x + Ur0 * gamma[n] / mu[n] * beta_sum

        #beta_check = np.exp(-mu[0] * (xi_diff))

        #pdb.set_trace()

        # eq. 1
        Cav = cn.sum(axis=1)

    else:

        cn = np.zeros(n_eigenmodes)

        for n in xrange(n_eigenmodes):

            # eq. 8
            beta[n, :] = np.exp(-mu[n] * (xi[-1] - xi[:]))

            # right hand side of eq. 7:
            beta_sum = np.sum((beta[n, 1:] - beta[n, :-1]) * Fa)

            # eq. 7
            cn[n] = x + Ur0 * gamma[n] / mu[n] * beta_sum

        # eq. 1
        Cav = cn.sum()

    # find age, not sure if Ur0 is the correct production term here...
    t_c = Cav / Ur0

    #
    #if t_c < decay_time / 2.0:
    #    print 'calculated age much smaller than decay time'
    #
    #    print 'need to solve age eq. numerically...., not implemented yet'

    return t_c


def calculate_he_age_meesters_dunai_2002(t, T, radius, U238, Th232,
                                         D0_div_a2=10**7.7,
                                         Ea=117.0e3,
                                         R=8.3144621,
                                         decay_constant_238U=4.916e-18,
                                         decay_constant_232Th=1.57e-18):

    """

    :param t:
    :param T:
    :param radius:
    :param U238:
    :param Th232:
    :param D0_div_a2:
    :param Ea:
    :param R:
    :param decay_constant_238U:
    :param decay_constant_232Th:
    :return:
    """

    D0 = D0_div_a2 * radius ** 2

    Dw = (D0 / radius**2 * np.exp(-Ea / (R*T))) * radius**2

    # diffusivity acc. to Cerniak et al (2009)
    Dc = 2.10e-6 * np.exp(-Ea / (R*T))

    # calculate He production rate
    Ur0 = 8 * U238 * decay_constant_238U + 6 * Th232 * decay_constant_232Th

    He_age = He_diffusion_Meesters_and_Dunai_2002(t, Dc, radius, Ur0, n_eigenmodes=15)

    return He_age