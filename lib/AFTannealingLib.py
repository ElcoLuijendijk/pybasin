"""
python library for apatite fission track annealing algorithms

Elco Luijendijk, Goettingen University

2011-2014

elco.luijendijk at geo.uni-goettingen.de

"""


import math
import sys
import pdb
import itertools
import numpy as np
from pylab import normpdf

# import fortran module 
try:
    import calculate_reduced_AFT_lengths
except ImportError:
    print('failed to import fortran annealing module')
    print('use slower python implementation of AFT annealing module instead')
    print('compile the fortran module by running the following command ' \
          'in the source directory of this module:')
    print('f2py -c calculate_reduced_AFT_lengths.f90 ' \
          '-m calculate_reduced_AFT_lengths')


def Cl_wt_fraction_to_APFU(Cl_wtfract):
    
    """
    calculate Cl atoms per formula unit (apfu)
    """
    
    Ca = 40.078
    P = 30.973762
    O = 15.9994
    H = 1.00794
    F = 18.9984032
    Cl = 35.453
    Br = 79.904

    #Cl_apfu = 0.33
    #F_apfu = 0.33
    #OH_apfu = 0.33

    # weight fluorapatite:
    fluorapatite_weight = Ca*10 + (P + 4*O)*6 + F*2
    
    Cl_apfu = fluorapatite_weight * Cl_wtfract / Cl
    
    return Cl_apfu


def convertDparFrom55to50M(Dpar50):
    """
    convert Dpar from Carlson et al (1999) 5.0 M etching conditions to
    Barbarand et al (2003) 5.5 M etching conditions
    
    source: Ketcham et al. (2007) Am. Min. 92: 799-810
    
    
    """
    
    Dpar55 = 0.9231 * Dpar50+0.2515
    
    return Dpar55


def calculate_normalized_density(r, r_crit=0.765, rmin=0.5275,
                                 d1=1.600, d2=-0.600, d3=9.205,
                                 d4=-9.157, d5=2.269):
    
    """
    convert reduced length (r) to normalized fission track density (rho)
    (eq. 7a and 7b in Ketcham 2000)
    
    input
    -----
        r       reduced fission track length
    
    optional keywords
    -----------------
        r_crit = 0.765
        d1 = 1.600
        d2 = -0.600
        d3 = 9.205
        d4 = -9.157
        d5 = 2.269
        rmin=0.5275
                cutoff value for track lengths. based on value reported by  
                Green (1988) Earth and Planetary Science Letters 89 (3-4)
        
    returns
    -------
        rho     normalized track density
    
    - checked nov 2012, equations are correct.
    - Not so sure about length cutoff limit, could be different than value
      used in HeFTy?
    
    """
    rho = np.zeros(r.shape)
    
    rho[r >= r_crit] = d1 * r[r >= r_crit] + d2
    rho[r < r_crit] = d3 * r[r < r_crit] ** 2 + d4 * r[r < r_crit] + d5
    
    # remove minimum r in Green's 1988 dataset of r vs track density
    rho[r <= rmin] = 0
    # remove <0 density values 
    rho[rho < 0] = 0 
    
    return rho


def caxis_project_reduced_lengths(rc, p1=-1.499, p2=4.150, p3=-1.656):

    """
    convert reduced lengths to c-axis projected reduced lengths
    see fig 8 in ketcham et al.,  1999
    
    Parameters
    ----------
    rc
        reduced fisssion track length, c-axis projected
    p1
        default = -1.499
    p2
        default = 4.150
    p3
        default = -1.656

    Returns
    -------
    rm
        non c-axis projected reduced length
    
    - checked eq., nov 2012, correct.
    
    """
    
    return p1 * (rc ** 2) + p2 * rc + p3


def get_initial_track_length(kinetic_parameter, kinetic_value,
                             apply_c_axis_correction, method='Carlson1999',
                             verbose=False):

    """
    set initial track lengths, 
    
    using linear correlations with kinetic parameters Dpar and Cl wt%
    following Carlson et al. (1999). Am. Min. 84, 1213-1223
    or Ketcham2007
    
    input
    -----
        kinetic_parameter        'Dpar' or 'Clwt'
        kinetic_value            value of kinetic parameter
        apply_c_axis_correction    option to calculate c-axis corrected i
                                   nitial track lengths
        method                  .......
        
    returns
    -------
        initial_track_length      initial fission track length (um)
        
    """
    
    if method == 'Carlson1999':
        if kinetic_parameter == 'Dpar':
            Dpar = kinetic_value
            if apply_c_axis_correction == False:
                initial_track_length = 15.63 + 0.283 * Dpar
            else:
                initial_track_length = 16.10 + 0.205 * Dpar
            
        elif kinetic_parameter == 'Clwt':
            Clwtfract = kinetic_value
            Cl_apfu = Cl_wt_fraction_to_APFU(Clwtfract)
            if apply_c_axis_correction == False:
                initial_track_length = 16.18 + 0.544 * Cl_apfu
                
            else:
                initial_track_length = 16.495 + 0.407 * Cl_apfu

    elif method == 'Ketcham2007':
        
        if kinetic_parameter == 'Dpar':
            Dpar = kinetic_value
            if apply_c_axis_correction  ==  False:
                initial_track_length = 15.391 + 0.258 * Dpar
            else:
                initial_track_length = 15.582 + 0.287 * Dpar
            
        elif kinetic_parameter == 'Clwt':
            Clwtfract = kinetic_value
            Cl_apfu = Cl_wt_fraction_to_APFU(Clwtfract)
            if apply_c_axis_correction is False:
                initial_track_length = 15.936 + 0.538 * Cl_apfu
                # use HeFTy 1.6.7 default values:
                initial_track_length2 = 15.936 + 0.16 * Clwtfract*100.0
            else:
                initial_track_length = 16.187 + 0.604 * Cl_apfu
                 # use HeFTy 1.6.7 default values:
                initial_track_length2 = 16.187 + 0.18 * Clwtfract*100.0

            if verbose is True:
                print('init. lengths: %0.3f, %0.3f' % (initial_track_length,
                                                       initial_track_length2))
        
        elif kinetic_parameter == 'rmr0':
            rmr0 = kinetic_value
            if apply_c_axis_correction is False:
                initial_track_length = 16.187
            else:
                initial_track_length = 15.936
            
    return initial_track_length


def calculate_kinetic_parameters(kinetic_parameter, kinetic_value):
    
    """
    set kinetic fission track annealling parameters
    according to Ketcham (2007)

    input
    -----
    
        kinetic_parameter        'Dpar' or 'Clwt'
        kinetic_value            value of kinetic parameter
    
    returns
    -------
        rmr0                    ...
        kappa                   ...
        
    """

    if kinetic_parameter == 'Clwt':
        # convert Chloride weight percentage to apfu:
        Cl_apfu = Cl_wt_fraction_to_APFU(kinetic_value)
        Clx = np.abs(Cl_apfu - 1.0)
        rmr0 = 0.83 * (((Clx - 0.13) / 0.87) ** 0.23)

    # calculate kinetic parameters (Dpar in um,  Cl in apfu):
    elif kinetic_parameter == 'Dpar':
        Dpar = kinetic_value
        rmr0 = 0.84 * (((4.58 - Dpar) / 2.98) ** 0.21)
        
    kappa = 1.04 - rmr0
    
    return rmr0, kappa



def set_annealing_parameters():
    """
    Ketcham et al. 2007 annealing model parameters
    
    returns:
        C0, C1, C2, C3, alpha, beta
    
    references:
        Ketcham et al. (2007) Am. Min. 92: 799-810,  Table 5c
    """
    
    # g function parameters:
    alpha = 0.04672
    beta = -1.0

    # fanning curvelinear function
    C0 = 0.39528
    C1 = 0.01073
    C2 = -65.12969
    C3 = -7.91715
    
    return C0, C1, C2, C3, alpha, beta


def correct_for_uranium_decay(time_bp, decay_const=1.551e-4):
    
    """
    correct for higher uranium conc in past 
    eq. 12 in Ketcham (2005), but modifed
    result checked for quoted 1% higher production at 64.2 Ma
    
    returns w: 
        weighting factor to correct for uranium decay over time
    inputs: 
        time: start en end of track formation,  in seconds bp:
    
    decay_const = 1.551e-4 
    total decay constant of ^238 U in Ma^-1
        =value from thesis Bart Hendriks.
    
    """

    w = np.exp(decay_const * time_bp[:-1])

    return w


def calculate_teq(g1, T1, C0, C1, C2, C3):
    
    """
    calculate equivalent time
    """
                  
    return math.exp(((g1 - C0) / C1 * (math.log(1.0 / T1) - C3)) + C2)
    
    
    
def calculate_reduced_track_lengths(dts, temperatures,
                                    alpha=0.04672,
                                    C0=0.39528, C1=0.01073,
                                    C2=-65.12969, C3=-7.91715):
    
    """
    Ketcham 2007 (HeFTy) annealing model
    
    input parameters:
        dts                     array of duration of each timestep (sec)
        temperatures            array of temperature during each timestep(K)
        alfa = 0.04672          emperically fitted annealing parameter,
                                see Ketcham et al. (2007) Am. Min.
        beta = -1.0             annealing parameter, see Ketcham(2007)
        C0 = 0.39528            annealing parameter, see Ketcham(2007)
        C1 = 0.01073            annealing parameter, see Ketcham(2007)
        C2 = -65.12969          annealing parameter, see Ketcham(2007)
        C3 = -7.91715           annealing parameter, see Ketcham(2007)
    
    returns:
        rc                  modeled reduced c-axis projected track length
    """

    #print '-' * 5
    #print 'start calculation of reduced track lengths'

    # initialize arrays:
    nsteps = len(dts)
    g = np.zeros(nsteps)
    rc = np.zeros(nsteps)
    dteq = np.zeros(nsteps)
    
    ##################################################################
    # calculate reduced track lengths
    ##################################################################  
    for j in range(nsteps):
        
        # reset g and dteq arrays
        g[:] = 0
        dteq[:] = 0
        
        # go through all timesteps and calculate equivalent time and g
        for ic, dt, T in zip(itertools.count(), dts[j:], temperatures[j:]):
            i = ic + j
            # calculate equivalent time (dteq)
            # equivalent time = time needed to reach annealing state of
            # previous timestep at current temperatures
            if i == j:
                dteq_i = 0
                dteq[i] = 0
            else:
                #h = ( (g[i-1] - C0) / C1 * (math.log(1.0/T)-C3)) + C2
                # dt[i-1] = e**h
                # h = ( (g[i-1] - C0) ) / (C1 * ln(1/T[i] - C3) )+ C2
                #dteq_i = calculate_teq(g[i-1], T, C0, C1, C2, C3)

                #math.exp(((g1 - C0) / C1 * (math.log(1.0 / T1) - C3)) + C2)

                dteq[i] = calculate_teq(g[i-1], T, C0, C1, C2, C3)

            # calculate g at each timestep:              
            g[i] = (C0 + C1 * ((np.log(dt + dteq[i]) - C2) / (np.log(1.0 / T)-C3)))

        #g = (C0 + C1 * ((np.log(dt + dteq) - C2) /
        #                    (np.log(1.0 / temperatures[j:])-C3)))

        # calculate reduced track length (rc) of tracks formed at timestep j
        # from value of g at the last timestep
        f = g[-1]**(1./alpha)
        rc[j] = 1.0 / (f + 1.0)
        #rc[j] = (1.0 - ((g[-1] * alpha +1.0)**(1.0/alpha)) * beta) **(1.0/beta)
        
    #pdb.set_trace()
        
    return rc


def kinetic_modifier_reduced_lengths(rc, rmr0, kappa):

    """
    correct modeled reduced track lengths for annealing resistance

    input parameters:
    rc                  c-axis projected reduced track length
    rmr0                annealing kinetics parameter
    kappa               annealing kinetics parameter

    returns:
    rc_corrected        corrected reduced track length
    """

    rc_mod = ((rc - rmr0) / (1.0-rmr0)) ** kappa

    return rc_mod


def kinetic_modifier_reduced_lengths_inverse(rc_mod, rmr0, kappa):

    """
    correct modeled reduced track lengths for annealing resistance

    input parameters:
    rc                  c-axis projected reduced track length
    rmr0                annealing kinetics parameter
    kappa               annealing kinetics parameter

    returns:
    rc_corrected        corrected reduced track length
    """

    rc = rc_mod ** (1.0/kappa) * (1.0 - rmr0) + rmr0

    return rc



def resample_time_temp_input(timesteps, temperature, max_temp_change=3.5):

    """
    check if no >3.5 degrees temperature change per timestep, 
    and resample using linear interpolation if yes
    
    Parameters
    ----------
    timesteps
        1D array containing time or duration of timesteps (My)
    temperature
        1D array of temperature (degr. C)
    max_temp_change
        maximum allowed change in temperature, default value=3.5
    
    Returns
    -------
    time_new
    temperature_new

    """
    
    time_new = timesteps.copy()
    temperature_new = temperature.copy()
    
    # check if temperature changes exceed limit:
    delta_T = temperature_new[1:] - temperature_new[:-1]
    counter = 0
    while (np.abs(delta_T)).max() > max_temp_change:
        
        # find locations where T change > 3.5 degrees:
        exc_loc = np.argmax(np.abs(delta_T)) + 1
        
        T_insert = (temperature_new[exc_loc-1]
                    + temperature_new[exc_loc]) / 2.0
        t_insert = (time_new[exc_loc-1] + time_new[exc_loc]) / 2.0
        temperature_new = np.insert(temperature_new, exc_loc, T_insert)
        time_new = np.insert(time_new, exc_loc, t_insert)
        
        # check new temperature array:
        delta_T = temperature_new[1:] - temperature_new[:-1]
        
        counter += 1
        if counter > 1000:
            print('error in temperature resampling function')
            print('changes in temperature too high')
            pdb.set_trace()
    
    return time_new, temperature_new


def simulate_AFT_annealing(timesteps, temperature_input, kinetic_value,
                           verbose=False,
                           method='Ketcham2007',
                           apply_c_axis_correction=False,
                           kinetic_parameter='Clwt',
                           initial_track_length=-99999,
                           binsize=0.25,
                           rmr0_min=0,
                           rmr0_max=0.85,
                           kappa=None,
                           min_length=2.18,
                           surpress_resampling=False,
                           use_fortran_algorithm=True,
                           annealing_eq='FC',
                           alpha=0.04672,
                           C0=0.39528,
                           C1=0.01073,
                           C2=-65.12969,
                           C3=-7.91715):
    
    
    """
    Forward modeling of apatite fission track ages and 
    track length distributions
    
    based on algorithms published by Ketcham et al (2000) and
    Ketcham (2005, 2007), see references below.
    
    Parameters
    ----------
    timesteps:
        1D array of time steps  (My)
    temperature_input:
        1D array of temperature (degr. C)
    kinetic_value:
        value of the kinetic parameter, the kinetic parameter is
        defined by optional argument *kinetic_parameter, see below
    method='Ketcham2007':
        'Ketcham2007' for Ketcham 2007 annealing algorithm 
        (default in HeFTy), 'Ketcham2000' for Ketcham (2000)
        algorithm (default in AFTsolve)
    kinetic_parameter='Clwt':
        parameter to use for calculation of initial track length
        and annealing kinetics.
        'Clwt': use the chloride weight fraction of apatite,
        'Dpar': use the Dpar parameter as a measure of the
                annealing properties of apatite
        'rmr0': use rmr0
    initial_track_length = -99999:
        use this parameter to speicify an initial track length
        the default value of -99999 the initial track length is
        based on the kinetic parameter, using equations by
        Ketcham (2007)
    binsize=0.25:
        Size of the track length bins for the modeled 
        probability density function (pdf) of track lengths, in um.
        Default value is 0.25 um. The pdf is calculated over a 
        range of 0 to 20 um.
    rmr0_min=0:
        minimum allowed value of the kinetic paramter rmr0. 
    rmr0_max = 0.85
        maximum allowed value of the kinetic parameter rmr0.
        The default value of 0.85 is the highest value of rmr0 
        reported by Carlson (1999)
    min_length=2.18:
        minimum value of track legnth (um). All lengths velow this value
        are disregarded.
        Value taken from ref.....
    surpress_resampling=False:
        option to surpress resampling of time and temperature arrays to avoid
        large changes in temperature
    
    use_fortran_algorithm=True:
        Use a separate Fortran function to calculate reduced
        fission track lengths.
        results in much shorter runtimes.

    Returns
    -------
    track_length_pdf : numpy array
        probability density function of track lengths
    AFTage_corrected_My:
        apatite fission track age (My)
    l_mean:
        mean track length (um)
    l_mean_std:
        standard deviation of mean track length (um)
    rm:
        reduced track length (um)
    rc:
        c-axis corrected reduced track length (um)
    rho_age:
        
    dt:
        time steps durations (sec)


    References
    ----------
    
    Ketcham, R. A., R. A. Donelick, and M. B. Donelick. 2000. 
        AFTSolve: A program for multi-kinetic modeling of apatite fission-track data.
        Geological Materials Research 2, no. 1 (May): 1-32.
        http://ammin.geoscienceworld.org/cgi/content/abstract/88/5-6/929.
        
    Ketcham, Richard A. 2005. 
        Forward and Inverse Modeling of Low-Temperature Thermochronometry Data.
        Reviews in Mineralogy and Geochemistry 58, no. 1: 275314.
        doi:10.2138/rmg.2005.58.11.
    
    Ketcham, R. A., A. Carter, R. A. Donelick, J. Barbarand, and A. J. Hurford. 2007.
        Improved modeling of fission-track annealing in apatite.
        American Mineralogist 92, no. 5-6: 799810.
        doi:10.2138/am.2007.2281.
    
    """

    Myr = (1.0e6 * 365.0 * 24.0 * 60.0 * 60.0)
    
    ####################################################################
    if verbose is True:
        print('-' * 20)
        print('T-t path:')
        print('duration = %0.1f My; mean,  min,  max T = %.0f, %.0f, %.0f'\
            % (timesteps.max(), temperature_input.mean(),
               temperature_input.min(), temperature_input.max()))
    
    # convert temperature units from degr. C to Kelvin:
    temperature = temperature_input + 273.15

    ##########################################################
    # check if no >3.5 degrees temperature change per timestep
    ##########################################################
    if verbose is True:
        print('resampling time steps')

    if surpress_resampling is False:
        timesteps, temperature = resample_time_temp_input(timesteps,
                                                          temperature)
    
    delta_T = temperature[1:] - temperature[:-1]
    if (abs(delta_T)).max() > 3.5 and surpress_resampling == False:
        max_loc = np.argmax(abs(delta_T))
        
        msg = 'warning %0.2f degr. change in T in timestep %s of %s My' % ((abs(delta_T)).max(),
                                                                           max_loc,
                                                                           timesteps[max_loc])
        raise ValueError(msg)
    
    #####################################
    # set AFT annealing model parameters:
    #####################################
    nsteps = len(temperature)
    r_standard = 0.893

    ###########################################################
    # calculate evolution of a fissions track over time,  and 
    # record final (last timestep) reduced length (r)
    # and c-axis non-corrected length (r_cmod):
    ###########################################################
    if verbose is True:
        print('calculating reduced track lengths')
    
    # get duration of each timestep in seconds
    dts = (timesteps[1:] - timesteps[:-1]) * 1.0e6 * 365 * 24 * 60 * 60
    nsteps = len(dts)
    
    # take midpoint values of temperature array:
    if verbose is True:
        print('taking midpoint values of temperature input array')
    temperature = (temperature[1:] + temperature[:-1]) / 2.0
    
    # get annealing kinetics:
    if kinetic_parameter != 'rmr0':
        rmr0, kappa = \
            calculate_kinetic_parameters(kinetic_parameter, kinetic_value)
    else:
        if verbose is True:
            print('using rmr0 as kinetic parameter')
        rmr0 = kinetic_value
        if kappa is None:
            kappa = 1.04 - rmr0

    if verbose is True:
        print('rmr0 = %0.3f, kappa = %0.3f' % (rmr0, kappa))

    if np.isnan(rmr0) is True or rmr0 <= rmr0_min:
        print('!! warning, rmr0 lower than minimum')
        print('!! %s = %0.3f' % (kinetic_parameter, kinetic_value))
        print('!! setting rmr0 to %0.3f' % rmr0_min)
        rmr0 = rmr0_min
        kappa = 1.04 - rmr0
    elif rmr0 > rmr0_max:
        print('!! warning, rmr0 value exceeds most resistant apatite in ' \
              'Carlson (1999) dataset')
        print('!! adjusting rmr0 from %0.3f to %0.3f' % (rmr0, rmr0_max))
        rmr0 = rmr0_max
        kappa = 1.04 - rmr0
        
    if use_fortran_algorithm is True:

        # fortran module for reduced track lengths:
        # call fortran module to calculate reduced fission track lengths
        if annealing_eq is 'FA':
            annealing_eq_f90 = 1
        elif annealing_eq is 'FC':
            annealing_eq_f90 = 2

        try:

            rcf = calculate_reduced_AFT_lengths.reduced_ln(
                dts, temperature,
                rmr0, kappa,
                annealing_eq_f90,
                alpha, C0, C1, C2, C3,
                nsteps)

            # rmf, rcf = calculate_reduced_AFT_lengths.reduced_ln(
            #    dts, temperature, rmr0, kappa, nsteps)
            rmf = caxis_project_reduced_lengths(rcf)

            # correct 0 length tracks:
            rmf[rmf < 0] = 0.0

            rm = rmf
            rc = rcf

        except NameError:
            if verbose is True:
                print('use python reduced track length function instead of fortran')
            r_cmod = calculate_reduced_track_lengths(dts, temperature,
                                                     C0=C0, C1=C1, C2=C2, C3=C3,
                                                     alpha=alpha)
            rcp = kinetic_modifier_reduced_lengths(r_cmod, rmr0, kappa)
            rmp = caxis_project_reduced_lengths(rcp)
            # correct 0 length tracks:
            rmp[rmp < 0] = 0.0
            rm = rmp
            rc = rcp

            #pdb.set_trace()

    else:
        if verbose is True:
            print('use python reduced track length function instead of fortran')
        # python reduced track length function:
        r_cmod = calculate_reduced_track_lengths(dts, temperature,
                                                 C0=C0, C1=C1, C2=C2, C3=C3,
                                                 alpha=alpha)
        rcp = kinetic_modifier_reduced_lengths(r_cmod, rmr0, kappa)
        rmp = caxis_project_reduced_lengths(rcp)
        # correct 0 length tracks:
        rmp[rmp < 0] = 0.0
        rm = rmp
        rc = rcp

    if verbose is True:
        print('final reduced lengths rm = %0.3f, rc = %0.3f' % (rm[-1], rc[-1]))

    #pdb.set_trace()

    ##########################################################
    # calculate weighting factor to correct for uranium decay 
    # (eq. 6 in Ketcham,  2000)
    ##########################################################
    time_Ma = timesteps.max() - timesteps
    w = correct_for_uranium_decay(time_Ma)    
    
    if verbose is True:
        print('w mean, min, max %0.2e, %0.2e, %0.2e'\
            % (w.mean(), w.min(), w.max()))
        
    ###################################################
    # calculate observation frequency:
    ###################################################
    rho = calculate_normalized_density(rc)
    
    if verbose is True:
        print('calculated observation frequency')
        print('rho mean, min, max %0.2e, %0.2e, %0.2e'\
            % (rho.mean(), rho.min(), rho.max()))
    
    ###########################################
    # set initial track lengths
    ##########################################
    l0 = initial_track_length
    if initial_track_length <= 0 or initial_track_length > 20:
        l0 = get_initial_track_length(kinetic_parameter,
                                      kinetic_value,
                                      apply_c_axis_correction,
                                      method=method)
    if verbose is True:
        print('initial track length %0.2f' % l0)
    
    #################################
    # calculate track lengths
    #################################
    l = rm * l0
    # remove short track lengths (<2.18 um)
    l = np.where(l < min_length, 0, 1) * l
   
    # calculate track length standard deviation,
    # using eqs. from fig 1 in Ketcham (2000)
    if apply_c_axis_correction is True:
        l_std = 0.008452 * l ** 2 - 0.2442 * l + 2.312
        l_std[np.where(l_std > 3.)] = 3.
    else:
        l_std = 0.02858 * l ** 2 - 0.8733 * l + 7.464
        l_std[np.where(l_std > 4.)] = 4.
        
    #######################################    
    # calculate probability density function
    #######################################
    if verbose is True:
        print('start calculation of pdf of track lengths:')
    track_ln_prob = np.zeros((nsteps, int(20 / binsize)))
    bins_ = np.arange(0, 20, binsize)
    for i in range(nsteps):
        track_ln_prob[i, :] = normpdf(bins_, l[i], l_std[i])
    
    #####################################################
    # correct for uranium decay and observation frequency
    #####################################################
    for j in range(0, int(20 / binsize)):
        track_ln_prob[:, j] = track_ln_prob[:, j] * w * rho
    
    # sum probability density,  and normalize
    if verbose is True:
        print('sum probability density,  and normalize')
    track_length_pdf = np.zeros((int(20 / binsize)))
    for j in range(0, int(20 / binsize)):
        track_length_pdf[j] = track_ln_prob[:, j].sum()
    track_length_pdf = track_length_pdf/track_length_pdf[:].sum()
    if verbose is True:
        print('done calculating probability density track lengths')
        
    ####################################################################
    # equation 15+16,  mean and standard deviation of model track length
    ####################################################################
    if verbose is True:
        print('calculate mean and std of track length from PDF')
    l_mean = 0
    for j in range(0, int(20 / binsize)):
        l_mean += (((j * binsize)+(0.5 * binsize)) * track_length_pdf[j])
    dummy = 0
    for j in range(len(track_length_pdf)):
        dummy += track_length_pdf[j] * ((j * binsize)-l_mean) ** 2
    l_mean_std = np.sqrt(dummy)
    l_median = (track_length_pdf[:].argmax() * binsize)

    if verbose is True:
        print('mean track length  =  %0.2f,  std =  %0.2f,  median: %0.3f'\
            % (l_mean, l_mean_std, l_median))
    if np.isnan(l_mean) is True:
        print('warning, track length calculation failed')
        if verbose is True:
            pdb.set_trace()
    
    ###################################################
    # calculate observation frequency:
    ###################################################
    rho_s = r_standard
    dt = dts

    #####################
    # calculate AFT ages:
    #####################
    # take midpoint value of reduced track length (rc)
    rc_mid = rc.copy()
    rc_mid[0] = rc[0]
    rc_mid[1:] = (rc[1:] + rc[:-1]) * 0.5
    
    # calculate fission track age density:
    rho_age = calculate_normalized_density(rc_mid) * w
    
    if verbose is True:
        print('dt mean, min, max %0.2e, %0.2e, %0.2e' \
              % (dt.mean(), dt.min(), dt.max()))
        print('rho_age mean, min, max %0.2e, %0.2e, %0.2e'\
              % (rho_age.mean(), rho_age.min(), rho_age.max()))
    
    aft_age_uncorrected = 0
    for i in range(nsteps):
        aft_age_uncorrected += dt[i] * rho_age[i]

    aft_age_corrected = aft_age_uncorrected / rho_s

    aft_age_myr = aft_age_corrected / Myr

    if verbose is True:
        print('AFT age = %0.2f My, avg rho = %0.3f, rho standard = %s'\
            % (aft_age_myr, rho_age.mean(), rho_s))

    if aft_age_corrected == 0:
        track_length_pdf[:] = 0
        l_mean = 0
        l_mean_std = 0

    return track_length_pdf, aft_age_myr, l_mean, l_mean_std, rm, rc, rho_age, dt
