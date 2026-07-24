"""
python library for apatite fission track annealing algorithms

Elco Luijendijk, Goettingen University

2011-2014

elco.luijendijk at geo.uni-goettingen.de

"""


import math
import sys
import itertools
import logging
import numpy as np
#from pylab import normpdf
import scipy.stats
from numba import jit

logger = logging.getLogger(__name__)


# import fortran module
try:
    import calculate_reduced_AFT_lengths
    logger.info('import fortran annealing module ok')
except:
    try:
        import lib.calculate_reduced_AFT_lengths as calculate_reduced_AFT_lengths
        logger.info('import fortran annealing module ok')
    except ImportError as e:
        logger.debug(e)
        logger.info('fortran annealing module not found, using slower pure-python '
                    'implementation (run "f2py -c calculate_reduced_AFT_lengths.f90 '
                    '-m calculate_reduced_AFT_lengths" to speed this up)')

@jit(nopython=True)
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

@jit(nopython=True)
def convertDparFrom55to50M(Dpar50):
    """
    convert Dpar from Carlson et al (1999) 5.0 M etching conditions to
    Barbarand et al (2003) 5.5 M etching conditions
    
    source: Ketcham et al. (2007) Am. Min. 92: 799-810
    
    
    """
    
    Dpar55 = 0.9231 * Dpar50+0.2515
    
    return Dpar55


@jit(nopython=True)
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


@jit(nopython=True)
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


@jit(nopython=True)
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

            #if verbose is True:
            #    print('init. lengths: %0.3f, %0.3f' % (initial_track_length,
            #                                           initial_track_length2))
        
        elif kinetic_parameter == 'rmr0':
            rmr0 = kinetic_value
            if apply_c_axis_correction is False:
                initial_track_length = 16.187
            else:
                initial_track_length = 15.936

    elif method == 'Ketcham2000':

        # the Ketcham et al. (1999, 2000) model predates the c-axis
        # projection technique introduced in Ketcham (2007),  so there is
        # no separate uncorrected calibration,  the same formula is used
        # regardless of apply_c_axis_correction

        if kinetic_parameter == 'Dpar':
            Dpar = kinetic_value
            # value used in HeFTy for the Ketcham (1999) model
            initial_track_length = 0.35 * Dpar + 15.72

        # no kinetic model specific to Cl content or rmr0 is available
        # for the Ketcham (1999, 2000) case,  the fixed intercept of the
        # Dpar based formula above is reused here as an approximation
        elif kinetic_parameter == 'Clwt':
            initial_track_length = 15.72

        elif kinetic_parameter == 'rmr0':
            initial_track_length = 15.72

    return initial_track_length


@jit(nopython=True)
def calculate_kinetic_parameters(kinetic_parameter, kinetic_value,
                                 method='Ketcham2007'):

    """
    set kinetic fission track annealling parameters,  following either
    Ketcham (2007) or Ketcham et al. (1999, 2000)

    input
    -----

        kinetic_parameter        'Dpar' or 'Clwt'
        kinetic_value            value of kinetic parameter
        method                  'Ketcham2007' (default) or 'Ketcham2000'

    returns
    -------
        rmr0                    ...
        kappa                   ...

    """

    if method == 'Ketcham2000':

        if kinetic_parameter == 'Dpar':
            # Ketcham et al. (1999), eq. 9a
            Dpar = kinetic_value
            rmr0 = 1.0 - np.exp(0.647 * (Dpar - 1.75) - 1.834)

        # no kinetic model specific to Cl content is available for
        # the Ketcham (1999, 2000) case,  the Ketcham (2007) conversion
        # is reused here as an approximation
        elif kinetic_parameter == 'Clwt':
            Cl_apfu = Cl_wt_fraction_to_APFU(kinetic_value)
            Clx = np.abs(Cl_apfu - 1.0)
            rmr0 = 0.83 * (((Clx - 0.13) / 0.87) ** 0.23)

        # Ketcham et al. (1999), eq. 8:  kappa = 1 - rmr0
        kappa = 1.0 - rmr0

    else:

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


@jit(nopython=True)
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


@jit(nopython=True)
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


@jit(nopython=True)
def calculate_teq(g1, T1, C0, C1, C2, C3):
    
    """
    calculate equivalent time
    """
                  
    return math.exp(((g1 - C0) / C1 * (math.log(1.0 / T1) - C3)) + C2)
    
    
@jit(nopython=True)  
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
        
        counter = list(range(len(dts[j:])))

        # go through all timesteps and calculate equivalent time and g
        #for ic, dt, T in zip(itertools.count(), dts[j:], temperatures[j:]):
        for ic, dt, T in zip(counter, dts[j:], temperatures[j:]):
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

                #dteq[i] = calculate_teq(g[i-1], T, C0, C1, C2, C3)
                dteq[i] = math.exp(((g[i-1] - C0) / C1 * (math.log(1.0 / T) - C3)) + C2)

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


@jit(nopython=True)
def calculate_reduced_track_lengths_matrix(dts, temperatures,
                                           alpha=0.04672,
                                           C0=0.39528, C1=0.01073,
                                           C2=-65.12969, C3=-7.91715):

    """
    Same recurrence as calculate_reduced_track_lengths, but keeps the
    reduced length reached at every intermediate timestep i for a
    track formed at timestep j, instead of only the value reached at
    the last timestep.

    Only needed to support all_timesteps=True in
    simulate_AFT_annealing, memory cost is O(nsteps**2) instead of
    O(nsteps).

    input parameters:
        dts                     array of duration of each timestep (sec)
        temperatures            array of temperature during each timestep(K)
        alfa = 0.04672          emperically fitted annealing parameter,
                                see Ketcham et al. (2007) Am. Min.
        C0 = 0.39528            annealing parameter, see Ketcham(2007)
        C1 = 0.01073            annealing parameter, see Ketcham(2007)
        C2 = -65.12969          annealing parameter, see Ketcham(2007)
        C3 = -7.91715           annealing parameter, see Ketcham(2007)

    returns:
        rc_matrix   nsteps by nsteps array, rc_matrix[j, i] is the
                    reduced length of a track formed at timestep j,
                    evaluated using the sub history from j through i
                    (i >= j), nan for i < j
    """

    nsteps = len(dts)
    g = np.zeros(nsteps)
    dteq = np.zeros(nsteps)
    rc_matrix = np.zeros((nsteps, nsteps))
    rc_matrix[:] = np.nan

    for j in range(nsteps):

        g[:] = 0
        dteq[:] = 0

        for i in range(j, nsteps):

            dt = dts[i]
            T = temperatures[i]

            if i == j:
                dteq[i] = 0
            else:
                dteq[i] = math.exp(((g[i - 1] - C0) / C1 *
                                   (math.log(1.0 / T) - C3)) + C2)

            g[i] = (C0 + C1 * ((np.log(dt + dteq[i]) - C2) /
                              (np.log(1.0 / T) - C3)))

            f = g[i] ** (1.0 / alpha)
            rc_matrix[j, i] = 1.0 / (f + 1.0)

    return rc_matrix


@jit(nopython=True)
def calculate_reduced_track_lengths_ketcham1999(dts, temperatures,
                                                alpha=-0.12327,
                                                beta=-11.988,
                                                C0=-19.844, C1=0.38951,
                                                C2=-51.253, C3=-7.6423):

    """
    Ketcham et al. (1999) fanning curvilinear annealing model,  as
    implemented in the AFTSolve program of Ketcham et al. (2000)

    Unlike calculate_reduced_track_lengths,  which follows Ketcham et al.
    (2007) with the equation shape parameter beta fixed at -1,  this
    model keeps beta as a free,  separately fitted parameter. Because of
    this,  the reduced length reached at the end of the previous timestep
    has to be converted back through the general length transform g(r)
    before it can be used to find the equivalent time at a new
    temperature,  rather than carrying a single transformed value forward
    as calculate_reduced_track_lengths does.

    Two guards against invalid,  fractional powers of a negative number
    are needed,  following Ketcham (1999, 2005) and matching the
    behaviour of HeFTy: a track already annealed below a reduced length
    of 0.0007 is treated as fully annealed and stays at a reduced length
    of 0 for all following timesteps; a new reduced length is only
    computed from the general formula while the intermediate term
    (alpha * f + 1) stays at or above 0.00002,  otherwise the track is
    also treated as fully annealed at the current timestep.

    input parameters:
        dts                     array of duration of each timestep (sec)
        temperatures            array of temperature during each timestep (K)
        alpha                   emperically fitted annealing parameter,
                                see Ketcham et al. (1999) Am. Min.
        beta                    annealing parameter,  see Ketcham et al. (1999)
        C0 = -19.844            annealing parameter,  see Ketcham et al. (1999)
        C1 = 0.38951            annealing parameter,  see Ketcham et al. (1999)
        C2 = -51.253            annealing parameter,  see Ketcham et al. (1999)
        C3 = -7.6423            annealing parameter,  see Ketcham et al. (1999)

    returns:
        rc                  modeled reduced c-axis projected track length
    """

    nsteps = len(dts)
    rc = np.zeros(nsteps)

    for j in range(nsteps):

        r = 0.0

        for ic in range(len(dts[j:])):

            i = ic + j
            dt = dts[i]
            T = temperatures[i]

            if i == j:
                # track formed at this timestep,  no prior annealing state
                dteq = 0.0
            elif r < 0.0007:
                # already fully annealed,  stays there permanently
                r = 0.0
                continue
            else:
                g_prev = (((1.0 - r ** beta) / beta) ** alpha - 1.0) / alpha
                dteq = math.exp(((g_prev - C0) / C1 *
                                 (math.log(1.0 / T) - C3)) + C2)

            f = (C0 + C1 * ((math.log(dt + dteq) - C2) /
                           (math.log(1.0 / T) - C3)))

            inner_base = alpha * f + 1.0
            if inner_base < 0.00002:
                r = 0.0
            else:
                inner_root = inner_base ** (1.0 / alpha)
                outer_base = 1.0 - beta * inner_root
                r = outer_base ** (1.0 / beta)

        rc[j] = r

    return rc


@jit(nopython=True)
def calculate_reduced_track_lengths_ketcham1999_matrix(dts, temperatures,
                                                        alpha=-0.12327,
                                                        beta=-11.988,
                                                        C0=-19.844, C1=0.38951,
                                                        C2=-51.253, C3=-7.6423):

    """
    Same recurrence as calculate_reduced_track_lengths_ketcham1999, but
    keeps the reduced length reached at every intermediate timestep i
    for a track formed at timestep j, instead of only the value
    reached at the last timestep. See that function for the guards
    against invalid, fractional powers of a negative number.

    Only needed to support all_timesteps=True in
    simulate_AFT_annealing, memory cost is O(nsteps**2) instead of
    O(nsteps).

    returns:
        rc_matrix   nsteps by nsteps array, rc_matrix[j, i] is the
                    reduced length of a track formed at timestep j,
                    evaluated using the sub history from j through i
                    (i >= j), nan for i < j
    """

    nsteps = len(dts)
    rc_matrix = np.zeros((nsteps, nsteps))
    rc_matrix[:] = np.nan

    for j in range(nsteps):

        r = 0.0

        for i in range(j, nsteps):

            dt = dts[i]
            T = temperatures[i]

            if i == j:
                # track formed at this timestep, no prior annealing state
                dteq = 0.0
            elif r < 0.0007:
                # already fully annealed, stays there permanently
                r = 0.0
                rc_matrix[j, i] = r
                continue
            else:
                g_prev = (((1.0 - r ** beta) / beta) ** alpha - 1.0) / alpha
                dteq = math.exp(((g_prev - C0) / C1 *
                                 (math.log(1.0 / T) - C3)) + C2)

            f = (C0 + C1 * ((math.log(dt + dteq) - C2) /
                           (math.log(1.0 / T) - C3)))

            inner_base = alpha * f + 1.0
            if inner_base < 0.00002:
                r = 0.0
            else:
                inner_root = inner_base ** (1.0 / alpha)
                outer_base = 1.0 - beta * inner_root
                r = outer_base ** (1.0 / beta)

            rc_matrix[j, i] = r

    return rc_matrix


def calculate_reduced_track_lengths_vectorized(dts, temperatures,
                                               alpha=0.04672,
                                               C0=0.39528, C1=0.01073,
                                               C2=-65.12969, C3=-7.91715,
                                               return_matrix=False):
    """
    Vectorized equivalent of calculate_reduced_track_lengths.

    Produces numerically identical results to the original O(n²) scalar loop
    by processing all starting-track indices j simultaneously as a batch at
    each inner timestep, rather than sequentially. Each inner step is a NumPy
    vector operation over all active tracks instead of a Python scalar loop.

    The sequential recurrence (dteq[i] depends on g[i-1]) is preserved: we
    advance all tracks together through each timestep i, which maintains the
    same operation order and therefore identical floating-point results.

    Parameters
    ----------
    dts : array_like
        Duration of each timestep (seconds).
    temperatures : array_like
        Temperature during each timestep (Kelvin).
    alpha, C0, C1, C2, C3 : float
        Ketcham 2007 annealing parameters.
    return_matrix : bool
        if True, return the reduced length reached at every
        intermediate timestep i for a track formed at timestep j
        (nsteps by nsteps array, nan for i < j), instead of only the
        value reached at the last timestep. Needed to support
        all_timesteps=True in simulate_AFT_annealing_vectorized, the
        g matrix this is derived from is already computed by the loop
        below regardless of this option, so requesting it costs no
        extra loop iterations, only the final elementwise conversion
        of the full matrix instead of its last column.

    Returns
    -------
    rc : numpy array
        Modeled reduced c-axis projected track length for each timestep,
        or the full rc_matrix described above if return_matrix is True.
    """
    nsteps = len(dts)

    # Pre-compute log terms — same values used in original inner loop
    log_inv_T = np.log(1.0 / temperatures)      # shape (nsteps,)
    log_dts = np.log(dts)                        # shape (nsteps,)

    # g[j, i] = annealing state of track formed at step j, evaluated at step i
    # We accumulate in a 2D array but only the lower triangle is active:
    # track j exists for i >= j.
    g = np.zeros((nsteps, nsteps))
    rc = np.zeros(nsteps)

    # Advance all tracks through each timestep i together.
    # At timestep i, tracks j=0..i are active.
    for i in range(nsteps):
        # Indices of all tracks formed at or before timestep i
        # j = 0, 1, ..., i
        active = slice(0, i + 1)         # tracks j <= i

        if i == 0:
            # First step: dteq is zero for all tracks (they were just formed)
            dteq_active = np.zeros(i + 1)
        else:
            # g_prev[j] = g[j, i-1] for active tracks
            g_prev = g[active, i - 1]    # shape (i,) — track j=i has no prev step
            # For j == i the track is born this step; dteq = 0.
            # For j < i we compute dteq from g[j, i-1].
            dteq_active = np.zeros(i + 1)
            if i > 0:
                g_prev_old = g[:i, i - 1]   # shape (i,), tracks j=0..i-1
                dteq_active[:i] = np.exp(
                    ((g_prev_old - C0) / C1 * (log_inv_T[i] - C3)) + C2
                )
            # dteq_active[i] stays 0 (track born at step i)

        # Compute g[j, i] for all active j
        g[active, i] = C0 + C1 * ((np.log(dts[i] + dteq_active) - C2)
                                   / (log_inv_T[i] - C3))

    if return_matrix is True:
        f_matrix = g ** (1.0 / alpha)
        rc_matrix = 1.0 / (f_matrix + 1.0)
        invalid = np.tril(np.ones((nsteps, nsteps), dtype=bool), k=-1)
        rc_matrix[invalid] = np.nan
        return rc_matrix

    # rc[j] is derived from g[j, nsteps-1] — the final annealing state
    g_final = g[:, nsteps - 1]
    f = g_final ** (1.0 / alpha)
    rc = 1.0 / (f + 1.0)

    return rc


def calculate_reduced_track_lengths_ketcham1999_vectorized(dts, temperatures,
                                                            alpha=-0.12327,
                                                            beta=-11.988,
                                                            C0=-19.844, C1=0.38951,
                                                            C2=-51.253, C3=-7.6423,
                                                            return_matrix=False):
    """
    Vectorized equivalent of calculate_reduced_track_lengths_ketcham1999.

    Produces numerically identical results to the scalar version by
    processing all starting-track indices j simultaneously as a batch at
    each inner timestep, rather than sequentially. Unlike
    calculate_reduced_track_lengths_vectorized (the Ketcham 2007 case),
    the reduced length r itself has to be carried forward and converted
    back through the general length transform g(r) at every timestep,
    since beta is a free parameter here rather than fixed at -1.

    Parameters
    ----------
    dts : array_like
        Duration of each timestep (seconds).
    temperatures : array_like
        Temperature during each timestep (Kelvin).
    alpha, beta, C0, C1, C2, C3 : float
        Ketcham et al. (1999) fanning curvilinear annealing parameters.
    return_matrix : bool
        if True, return the reduced length reached at every
        intermediate timestep i for a track formed at timestep j
        (nsteps by nsteps array, nan for i < j), instead of only the
        value reached at the last timestep. Needed to support
        all_timesteps=True in simulate_AFT_annealing_vectorized, the r
        matrix this returns is already computed by the loop below
        regardless of this option, so requesting it costs no extra
        loop iterations.

    Returns
    -------
    rc : numpy array
        Modeled reduced c-axis projected track length for each timestep,
        or the full r matrix described above if return_matrix is True.
    """
    nsteps = len(dts)

    log_inv_T = np.log(1.0 / temperatures)

    # r[j, i] = reduced length of track formed at step j, evaluated at step i
    r = np.zeros((nsteps, nsteps))

    for i in range(nsteps):

        # tracks j = 0 .. i are active at this timestep
        n_active = i + 1

        if i == 0:
            # first step: only the freshly formed track j=0 exists,
            # no prior annealing state
            dteq_active = np.zeros(1)
        else:
            r_prev = r[:i, i - 1]

            # tracks already below the fully annealed threshold stay there
            not_annealed = r_prev >= 0.0007

            dteq_active = np.zeros(n_active)
            g_prev = np.zeros(i)
            g_prev[not_annealed] = (((1.0 - r_prev[not_annealed] ** beta) / beta)
                                    ** alpha - 1.0) / alpha
            dteq_active[:i][not_annealed] = np.exp(
                ((g_prev[not_annealed] - C0) / C1 * (log_inv_T[i] - C3)) + C2)
            # dteq_active[i] stays 0 (track formed at this timestep)
            # dteq_active[:i][~not_annealed] stays 0 too, those tracks are
            # forced back to a fully annealed state below regardless

        f = C0 + C1 * ((np.log(dts[i] + dteq_active) - C2) / (log_inv_T[i] - C3))

        inner_base = alpha * f + 1.0
        fully_annealed_now = inner_base < 0.00002

        inner_root = np.zeros(n_active)
        inner_root[~fully_annealed_now] = (
            inner_base[~fully_annealed_now] ** (1.0 / alpha))
        outer_base = 1.0 - beta * inner_root

        r_new = np.zeros(n_active)
        r_new[~fully_annealed_now] = outer_base[~fully_annealed_now] ** (1.0 / beta)

        if i > 0:
            # tracks already fully annealed coming into this step stay at
            # exactly zero, regardless of the computation above
            r_new[:i][~not_annealed] = 0.0

        r[:n_active, i] = r_new

    if return_matrix is True:
        invalid = np.tril(np.ones((nsteps, nsteps), dtype=bool), k=-1)
        r[invalid] = np.nan
        return r

    rc = r[:, nsteps - 1]

    return rc


@jit(nopython=True)
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

    base = (rc - rmr0) / (1.0 - rmr0)

    # tracks annealed below the kinetic resistance threshold rmr0 are
    # fully annealed, clip to zero to avoid raising a negative base to
    # a fractional power (kappa is not an integer in general)
    base = np.where(base < 0.0, 0.0, base)

    rc_mod = base ** kappa

    return rc_mod


@jit(nopython=True)
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


def calculate_age_evolution(rc_matrix, dts, timesteps, rmr0, kappa,
                            rho_s=0.893):

    """
    Reconstruct the apatite fission track age the model would report
    if the thermal history had stopped, and the sample had been
    collected, right after each timestep, from the full rc_matrix
    computed by calculate_reduced_track_lengths_matrix or
    calculate_reduced_track_lengths_ketcham1999_matrix.

    A plain running sum of the age density calculated over the full
    thermal history is not equivalent to this quantity, and gives
    wrong, in some cases zero, ages for thermal histories that involve
    reheating. This is because the reduced length of a track formed at
    timestep j depends on the entire sub history from j through the
    last timestep of the run, so a track's contribution to the age
    density changes depending on how far the history is followed.
    Recomputing the uranium decay weighting factor relative to each
    intermediate timestep, from the matching column of rc_matrix,
    gives the same age that a full rerun truncated at that timestep
    would produce.

    input parameters:
        rc_matrix       nsteps by nsteps array from
                        calculate_reduced_track_lengths_matrix or the
                        Ketcham1999 equivalent
        dts             array of duration of each timestep (sec)
        timesteps       array of time values (My), length nsteps + 1
        rmr0, kappa     kinetic parameters, see
                        calculate_kinetic_parameters
        rho_s = 0.893   standard track density

    returns:
        aft_ages_myr    array of length nsteps, aft_ages_myr[k - 1] is
                        the AFT age (My) using only the first k
                        timesteps of the thermal history,
                        aft_ages_myr[-1] equals the age from the full
                        history
    """

    Myr = (1.0e6 * 365.0 * 24.0 * 60.0 * 60.0)
    nsteps = len(dts)
    aft_ages_myr = np.zeros(nsteps)

    for k in range(1, nsteps + 1):

        rc = kinetic_modifier_reduced_lengths(rc_matrix[:k, k - 1], rmr0, kappa)
        rc = np.where(rc < 0, 0.0, rc)

        rc_mid = rc.copy()
        rc_mid[1:] = (rc[1:] + rc[:-1]) * 0.5

        time_bp = timesteps[k] - timesteps[:k + 1]
        w = correct_for_uranium_decay(time_bp)

        rho_age = calculate_normalized_density(rc_mid) * w

        aft_age_uncorrected = np.sum(dts[:k] * rho_age)
        aft_ages_myr[k - 1] = aft_age_uncorrected / rho_s / Myr

    return aft_ages_myr


#@jit(nopython=True)
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
            logger.error('error in temperature resampling function')
            logger.info('changes in temperature too high')
    
    return time_new, temperature_new


#@jit(nopython=True)
def simulate_AFT_annealing(timesteps, temperature_input, kinetic_value,
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
                           alpha=None,
                           beta=None,
                           C0=None,
                           C1=None,
                           C2=None,
                           C3=None,
                           all_timesteps=False,
                           verbose=False):


    """
    Forward modeling of apatite fission track ages and
    track length distributions

    based on algorithms published by Ketcham et al (1999, 2000) and
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
        'Ketcham2007' for the Ketcham et al. (2007) fanning
        curvilinear annealing algorithm (default in HeFTy), with the
        equation shape parameter beta fixed at -1.
        'Ketcham2000' for the Ketcham et al. (1999) fanning
        curvilinear algorithm as implemented in the AFTSolve program
        of Ketcham et al. (2000), with beta as a free, separately
        fitted parameter. This option always uses the pure Python
        annealing function, use_fortran_algorithm is ignored, since
        no compiled Fortran implementation of this model exists.
    alpha, beta, C0, C1, C2, C3 = None:
        annealing model parameters, defaults depend on method: the
        Ketcham (2007) values (alpha=0.04672, beta=-1, C0=0.39528,
        C1=0.01073, C2=-65.12969, C3=-7.91715) for method='Ketcham2007',
        or the Ketcham (1999) values (alpha=-0.12327, beta=-11.988,
        C0=-19.844, C1=0.38951, C2=-51.253, C3=-7.6423) for
        method='Ketcham2000'. Passing an explicit value overrides only
        that parameter, the rest keep their method-specific default.
    kinetic_parameter='Clwt':
        parameter to use for calculation of initial track length
        and annealing kinetics.
        'Clwt': use the chloride weight fraction of apatite,
        'Dpar': use the Dpar parameter as a measure of the
                annealing properties of apatite
        'rmr0': use rmr0
    apply_c_axis_correction=False:
        option to use c-axis corrected initial track lengths and
        length standard deviations for method='Ketcham2007'. Ignored
        for method='Ketcham2000': that model's reduced length is
        already c-axis projected, ages and mean lengths come out
        identical regardless of this option, only the length standard
        deviation formula still differs between True and False.
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
    all_timesteps=False:
        option to also return the AFT age evaluated at every
        intermediate timestep, i.e. the age the model would report if
        the thermal history had stopped, and the sample had been
        collected, at that point, rather than only the age for the
        full input history. When True, an additional array is
        appended to the return values, see aft_ages_myr below. Needs
        the pure Python annealing function: use
        use_fortran_algorithm=False, or method='Ketcham2000', which
        always uses the pure Python function. Raises ValueError if
        combined with use_fortran_algorithm=True and
        method='Ketcham2007'. The extra bookkeeping needed for this
        option is O(nsteps**2) in both runtime and memory, instead of
        the O(nsteps) cost of the default single final age, so leave
        this at its default value of False unless the age evolution is
        actually needed.

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
    aft_ages_myr:
        only returned if all_timesteps=True: array of AFT ages (My),
        one for every timestep, see all_timesteps above


    References
    ----------
    
    Ketcham, R. A., R. A. Donelick, and W. D. Carlson. 1999.
        Variability of apatite fission-track annealing kinetics: III.
        Extrapolation to geological time scales.
        American Mineralogist 84, no. 9: 1235-1255.
        doi:10.2138/am-1999-0903.

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

    # fill in method specific annealing model defaults for any
    # parameter the caller did not explicitly override
    if method == 'Ketcham2000':
        if alpha is None:
            alpha = -0.12327
        if beta is None:
            beta = -11.988
        if C0 is None:
            C0 = -19.844
        if C1 is None:
            C1 = 0.38951
        if C2 is None:
            C2 = -51.253
        if C3 is None:
            C3 = -7.6423
    else:
        if alpha is None:
            alpha = 0.04672
        if beta is None:
            beta = -1.0
        if C0 is None:
            C0 = 0.39528
        if C1 is None:
            C1 = 0.01073
        if C2 is None:
            C2 = -65.12969
        if C3 is None:
            C3 = -7.91715

    ####################################################################
    if verbose is True:
        logger.info('-' * 20)
        logger.info('T-t path:')
        logger.info('duration = %0.1f My; mean,  min,  max T = %.0f, %.0f, %.0f' % (timesteps.max(), temperature_input.mean(), temperature_input.min(), temperature_input.max()))

    # convert temperature units from degr. C to Kelvin:
    temperature = temperature_input + 273.15

    ##########################################################
    # check if no >3.5 degrees temperature change per timestep
    ##########################################################
    if surpress_resampling is False:

        if verbose is True:
            logger.info('resampling time steps')

        timesteps, temperature = resample_time_temp_input(timesteps,
                                                          temperature)

    delta_T = temperature[1:] - temperature[:-1]
    if (abs(delta_T)).max() > 3.5 and surpress_resampling == False:
        max_loc = np.argmax(abs(delta_T))
        
        msg = ('temperature change of %0.2f degrees C in a single timestep (at %s My) '
              'exceeds the 3.5 degree limit assumed by the fission track annealing model; '
              'resampling of the time-temperature path (surpress_resampling=False) should '
              'normally prevent this, check the input time-temperature history' % (
                  (abs(delta_T)).max(), timesteps[max_loc]))
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
        logger.info('calculating reduced track lengths')
    
    # get duration of each timestep in seconds
    dts = (timesteps[1:] - timesteps[:-1]) * Myr
    nsteps = len(dts)
    
    # take midpoint values of temperature array:
    if verbose is True:
        logger.info('taking midpoint values of temperature input array')
    temperature = (temperature[1:] + temperature[:-1]) / 2.0
    
    # kappa fallback constant is model specific: 1.0 - rmr0 for the
    # Ketcham (1999, 2000) model,  1.04 - rmr0 for Ketcham (2007)
    kappa_offset = 1.0 if method == 'Ketcham2000' else 1.04

    # get annealing kinetics:
    if kinetic_parameter != 'rmr0':
        rmr0, kappa = \
            calculate_kinetic_parameters(kinetic_parameter, kinetic_value,
                                         method=method)
    else:
        if verbose is True:
            logger.info('using rmr0 as kinetic parameter')
        rmr0 = kinetic_value
        if kappa is None:
            kappa = kappa_offset - rmr0

    if verbose is True:
        logger.info('rmr0 = %0.3f, kappa = %0.3f' % (rmr0, kappa))

    if np.isnan(rmr0) is True or rmr0 <= rmr0_min:
        logger.warning('!! warning, rmr0 lower than minimum')
        logger.info('!! %s = %0.3f' % (kinetic_parameter, kinetic_value))
        logger.info('!! setting rmr0 to %0.3f' % rmr0_min)
        rmr0 = rmr0_min
        kappa = kappa_offset - rmr0
    elif rmr0 > rmr0_max:
        logger.warning('!! warning, rmr0 value exceeds most resistant apatite in Carlson (1999) dataset')
        logger.info('!! adjusting rmr0 from %0.3f to %0.3f' % (rmr0, rmr0_max))
        rmr0 = rmr0_max
        kappa = kappa_offset - rmr0

    if method == 'Ketcham2000':

        # no compiled Fortran implementation of the free beta,
        # Ketcham (1999, 2000) model exists, use_fortran_algorithm is
        # ignored and the pure Python function is always used
        if verbose is True:
            logger.info('method=Ketcham2000: using pure Python annealing '
                       'function, use_fortran_algorithm is ignored')
        r_cmod = calculate_reduced_track_lengths_ketcham1999(
            dts, temperature, alpha=alpha, beta=beta,
            C0=C0, C1=C1, C2=C2, C3=C3)
        rcp = kinetic_modifier_reduced_lengths(r_cmod, rmr0, kappa)

        # unlike the Ketcham (2007) model, the Ketcham (1999, 2000)
        # reduced length as calibrated here is already the c-axis
        # projected quantity, caxis_project_reduced_lengths is not
        # applied
        rcp[rcp < 0] = 0.0
        rm = rcp
        rc = rcp

    elif use_fortran_algorithm is True:

        # fortran module for reduced track lengths:
        # call fortran module to calculate reduced fission track lengths
        if annealing_eq == 'FA':
            annealing_eq_f90 = 1
        elif annealing_eq == 'FC':
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
                logger.info('use python reduced track length function instead of fortran')
            r_cmod = calculate_reduced_track_lengths(dts, temperature,
                                                     C0=C0, C1=C1, C2=C2, C3=C3,
                                                     alpha=alpha)
            rcp = kinetic_modifier_reduced_lengths(r_cmod, rmr0, kappa)
            rmp = caxis_project_reduced_lengths(rcp)
            # correct 0 length tracks:
            rmp[rmp < 0] = 0.0
            rm = rmp
            rc = rcp

    else:
        if verbose is True:
            logger.info('use python reduced track length function instead of fortran')
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
        logger.info('final reduced lengths rm = %0.3f, rc = %0.3f' % (rm[-1], rc[-1]))

    ####################################################################
    # optionally also compute the AFT age at every intermediate timestep
    ####################################################################
    aft_ages_myr = None
    if all_timesteps is True:

        if method == 'Ketcham2000':
            rc_matrix = calculate_reduced_track_lengths_ketcham1999_matrix(
                dts, temperature, alpha=alpha, beta=beta,
                C0=C0, C1=C1, C2=C2, C3=C3)
        elif use_fortran_algorithm is True:
            msg = ('all_timesteps=True needs the reduced length reached at '
                   'every intermediate timestep, which the compiled '
                   'Fortran algorithm does not return; use '
                   'use_fortran_algorithm=False or method=\'Ketcham2000\'')
            raise ValueError(msg)
        else:
            rc_matrix = calculate_reduced_track_lengths_matrix(
                dts, temperature, alpha=alpha, C0=C0, C1=C1, C2=C2, C3=C3)

        aft_ages_myr = calculate_age_evolution(
            rc_matrix, dts, timesteps, rmr0, kappa)

    ##########################################################
    # calculate weighting factor to correct for uranium decay
    # (eq. 6 in Ketcham,  2000)
    ##########################################################
    time_Ma = timesteps.max() - timesteps
    w = correct_for_uranium_decay(time_Ma)    
    
    if verbose is True:
        logger.info('w mean, min, max %0.2e, %0.2e, %0.2e' % (w.mean(), w.min(), w.max()))
        
    ###################################################
    # calculate observation frequency:
    ###################################################
    rho = calculate_normalized_density(rc)
    
    if verbose is True:
        logger.info('calculated observation frequency')
        logger.info('rho mean, min, max %0.2e, %0.2e, %0.2e' % (rho.mean(), rho.min(), rho.max()))
    
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
        logger.info('initial track length %0.2f' % l0)
    
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
        logger.info('start calculation of pdf of track lengths:')
    track_ln_prob = np.zeros((nsteps, int(20 / binsize)))
    bins_ = np.arange(0, 20, binsize)
    for i in range(nsteps):
        #track_ln_prob[i, :] = normpdf(bins_, l[i], l_std[i])
        track_ln_prob[i, :] = scipy.stats.norm.pdf(bins_, l[i], l_std[i])

    #####################################################
    # correct for uranium decay and observation frequency
    #####################################################
    for j in range(0, int(20 / binsize)):
        track_ln_prob[:, j] = track_ln_prob[:, j] * w * rho
    
    # sum probability density,  and normalize
    if verbose is True:
        logger.info('sum probability density,  and normalize')
    track_length_pdf = np.zeros((int(20 / binsize)))
    for j in range(0, int(20 / binsize)):
        track_length_pdf[j] = track_ln_prob[:, j].sum()
    track_length_pdf = track_length_pdf/track_length_pdf[:].sum()
    if verbose is True:
        logger.info('done calculating probability density track lengths')
        
    ####################################################################
    # equation 15+16,  mean and standard deviation of model track length
    ####################################################################
    if verbose is True:
        logger.info('calculate mean and std of track length from PDF')
    l_mean = 0
    for j in range(0, int(20 / binsize)):
        l_mean += (((j * binsize)+(0.5 * binsize)) * track_length_pdf[j])
    dummy = 0
    for j in range(len(track_length_pdf)):
        dummy += track_length_pdf[j] * ((j * binsize)-l_mean) ** 2
    l_mean_std = np.sqrt(dummy)
    l_median = (track_length_pdf[:].argmax() * binsize)

    if verbose is True:
        logger.info('mean track length  =  %0.2f,  std =  %0.2f,  median: %0.3f' % (l_mean, l_mean_std, l_median))
    if np.isnan(l_mean) is True:
        logger.warning('track length calculation failed')
    
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
        logger.info('dt mean, min, max %0.2e, %0.2e, %0.2e' % (dt.mean() / Myr, dt.min() / Myr, dt.max() / Myr))
        logger.info('rho_age mean, min, max %0.2e, %0.2e, %0.2e' % (rho_age.mean(), rho_age.min(), rho_age.max()))
        logger.info('dt sum = %0.2f' % (dt.sum() / Myr))

    aft_age_uncorrected = 0
    for i in range(nsteps):
        aft_age_uncorrected += dt[i] * rho_age[i]

    aft_age_corrected = aft_age_uncorrected / rho_s

    aft_age_myr = aft_age_corrected / Myr

    if verbose is True:
        logger.info('AFT age = %0.2f My, avg rho = %0.3f, rho standard = %s' % (aft_age_myr, rho_age.mean(), rho_s))

    if aft_age_corrected == 0:
        track_length_pdf[:] = 0
        l_mean = 0
        l_mean_std = 0

    if all_timesteps is True:
        return (track_length_pdf, aft_age_myr, l_mean, l_mean_std, rm, rc,
                rho_age, dt, aft_ages_myr)

    return track_length_pdf, aft_age_myr, l_mean, l_mean_std, rm, rc, rho_age, dt


def simulate_AFT_annealing_vectorized(timesteps, temperature_input, kinetic_value,
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
                                      alpha=None,
                                      beta=None,
                                      C0=None,
                                      C1=None,
                                      C2=None,
                                      C3=None,
                                      all_timesteps=False,
                                      verbose=False):
    """
    Vectorized equivalent of simulate_AFT_annealing.

    Identical to simulate_AFT_annealing in all respects except:

    1. Uses calculate_reduced_track_lengths_vectorized (Item 1) instead of
       calculate_reduced_track_lengths when the Fortran module is unavailable.
       This replaces the O(n²) Python scalar loop with O(n) NumPy vector steps.

    2. The track-length PDF construction (Item 5) is fully vectorized:
       - norm.pdf computed over all timesteps × bins in one broadcast call
       - weight/density correction applied with array multiplication
       - bin summation and normalisation via np.sum(axis=0)
       - l_mean and l_mean_std computed with np.dot instead of scalar loops

    All other logic, parameters, and return values are identical to the
    original function, including all_timesteps: see
    simulate_AFT_annealing for details. The reduced length matrix
    needed for all_timesteps=True is computed with
    return_matrix=True on the vectorized helper functions, which costs
    no extra loop iterations there, unlike the scalar version where a
    separate matrix-returning function has to be called.
    """

    Myr = (1.0e6 * 365.0 * 24.0 * 60.0 * 60.0)

    # fill in method specific annealing model defaults for any
    # parameter the caller did not explicitly override
    if method == 'Ketcham2000':
        if alpha is None:
            alpha = -0.12327
        if beta is None:
            beta = -11.988
        if C0 is None:
            C0 = -19.844
        if C1 is None:
            C1 = 0.38951
        if C2 is None:
            C2 = -51.253
        if C3 is None:
            C3 = -7.6423
    else:
        if alpha is None:
            alpha = 0.04672
        if beta is None:
            beta = -1.0
        if C0 is None:
            C0 = 0.39528
        if C1 is None:
            C1 = 0.01073
        if C2 is None:
            C2 = -65.12969
        if C3 is None:
            C3 = -7.91715

    if verbose is True:
        logger.info('-' * 20)
        logger.info('T-t path:')
        logger.info('duration = %0.1f My; mean,  min,  max T = %.0f, %.0f, %.0f' % (timesteps.max(), temperature_input.mean(), temperature_input.min(), temperature_input.max()))

    # convert temperature units from degr. C to Kelvin:
    temperature = temperature_input + 273.15

    ##########################################################
    # check if no >3.5 degrees temperature change per timestep
    ##########################################################
    if surpress_resampling is False:
        if verbose is True:
            logger.info('resampling time steps')
        timesteps, temperature = resample_time_temp_input(timesteps, temperature)

    delta_T = temperature[1:] - temperature[:-1]
    if (abs(delta_T)).max() > 3.5 and surpress_resampling == False:
        max_loc = np.argmax(abs(delta_T))
        msg = ('temperature change of %0.2f degrees C in a single timestep (at %s My) '
              'exceeds the 3.5 degree limit assumed by the fission track annealing model; '
              'resampling of the time-temperature path (surpress_resampling=False) should '
              'normally prevent this, check the input time-temperature history' % (
                  (abs(delta_T)).max(), timesteps[max_loc]))
        raise ValueError(msg)

    #####################################
    # set AFT annealing model parameters:
    #####################################
    nsteps = len(temperature)
    r_standard = 0.893

    if verbose is True:
        logger.info('calculating reduced track lengths')

    # get duration of each timestep in seconds
    dts = (timesteps[1:] - timesteps[:-1]) * Myr
    nsteps = len(dts)

    # take midpoint values of temperature array:
    if verbose is True:
        logger.info('taking midpoint values of temperature input array')
    temperature = (temperature[1:] + temperature[:-1]) / 2.0

    # kappa fallback constant is model specific: 1.0 - rmr0 for the
    # Ketcham (1999, 2000) model,  1.04 - rmr0 for Ketcham (2007)
    kappa_offset = 1.0 if method == 'Ketcham2000' else 1.04

    # get annealing kinetics:
    if kinetic_parameter != 'rmr0':
        rmr0, kappa = \
            calculate_kinetic_parameters(kinetic_parameter, kinetic_value,
                                         method=method)
    else:
        if verbose is True:
            logger.info('using rmr0 as kinetic parameter')
        rmr0 = kinetic_value
        if kappa is None:
            kappa = kappa_offset - rmr0

    if verbose is True:
        logger.info('rmr0 = %0.3f, kappa = %0.3f' % (rmr0, kappa))

    if np.isnan(rmr0) is True or rmr0 <= rmr0_min:
        logger.warning('!! warning, rmr0 lower than minimum')
        logger.info('!! %s = %0.3f' % (kinetic_parameter, kinetic_value))
        logger.info('!! setting rmr0 to %0.3f' % rmr0_min)
        rmr0 = rmr0_min
        kappa = kappa_offset - rmr0
    elif rmr0 > rmr0_max:
        logger.warning('!! warning, rmr0 value exceeds most resistant apatite in Carlson (1999) dataset')
        logger.info('!! adjusting rmr0 from %0.3f to %0.3f' % (rmr0, rmr0_max))
        rmr0 = rmr0_max
        kappa = kappa_offset - rmr0

    if method == 'Ketcham2000':

        # no compiled Fortran implementation of the free beta,
        # Ketcham (1999, 2000) model exists, use_fortran_algorithm is
        # ignored and the vectorized pure Python function is always used
        if verbose is True:
            logger.info('method=Ketcham2000: using vectorized pure Python '
                       'annealing function, use_fortran_algorithm is ignored')
        r_cmod = calculate_reduced_track_lengths_ketcham1999_vectorized(
            dts, temperature, alpha=alpha, beta=beta,
            C0=C0, C1=C1, C2=C2, C3=C3)
        rcp = kinetic_modifier_reduced_lengths(r_cmod, rmr0, kappa)

        # unlike the Ketcham (2007) model, the Ketcham (1999, 2000)
        # reduced length as calibrated here is already the c-axis
        # projected quantity, caxis_project_reduced_lengths is not
        # applied
        rcp[rcp < 0] = 0.0
        rm = rcp
        rc = rcp

    elif use_fortran_algorithm is True:
        if annealing_eq == 'FA':
            annealing_eq_f90 = 1
        elif annealing_eq == 'FC':
            annealing_eq_f90 = 2

        try:
            rcf = calculate_reduced_AFT_lengths.reduced_ln(
                dts, temperature,
                rmr0, kappa,
                annealing_eq_f90,
                alpha, C0, C1, C2, C3,
                nsteps)
            rmf = caxis_project_reduced_lengths(rcf)
            rmf[rmf < 0] = 0.0
            rm = rmf
            rc = rcf

        except NameError:
            if verbose is True:
                logger.info('use vectorized python reduced track length function')
            # --- Item 1: use vectorized fallback instead of scalar O(n²) loop ---
            r_cmod = calculate_reduced_track_lengths_vectorized(
                dts, temperature, C0=C0, C1=C1, C2=C2, C3=C3, alpha=alpha)
            rcp = kinetic_modifier_reduced_lengths(r_cmod, rmr0, kappa)
            rmp = caxis_project_reduced_lengths(rcp)
            rmp[rmp < 0] = 0.0
            rm = rmp
            rc = rcp

    else:
        if verbose is True:
            logger.info('use vectorized python reduced track length function')
        # --- Item 1: use vectorized fallback instead of scalar O(n²) loop ---
        r_cmod = calculate_reduced_track_lengths_vectorized(
            dts, temperature, C0=C0, C1=C1, C2=C2, C3=C3, alpha=alpha)
        rcp = kinetic_modifier_reduced_lengths(r_cmod, rmr0, kappa)
        rmp = caxis_project_reduced_lengths(rcp)
        rmp[rmp < 0] = 0.0
        rm = rmp
        rc = rcp

    if verbose is True:
        logger.info('final reduced lengths rm = %0.3f, rc = %0.3f' % (rm[-1], rc[-1]))

    ####################################################################
    # optionally also compute the AFT age at every intermediate timestep
    ####################################################################
    aft_ages_myr = None
    if all_timesteps is True:

        if method == 'Ketcham2000':
            rc_matrix = calculate_reduced_track_lengths_ketcham1999_vectorized(
                dts, temperature, alpha=alpha, beta=beta,
                C0=C0, C1=C1, C2=C2, C3=C3, return_matrix=True)
        elif use_fortran_algorithm is True:
            msg = ('all_timesteps=True needs the reduced length reached at '
                   'every intermediate timestep, which the compiled '
                   'Fortran algorithm does not return; use '
                   'use_fortran_algorithm=False or method=\'Ketcham2000\'')
            raise ValueError(msg)
        else:
            rc_matrix = calculate_reduced_track_lengths_vectorized(
                dts, temperature, alpha=alpha, C0=C0, C1=C1, C2=C2, C3=C3,
                return_matrix=True)

        aft_ages_myr = calculate_age_evolution(
            rc_matrix, dts, timesteps, rmr0, kappa)

    ##########################################################
    # calculate weighting factor to correct for uranium decay
    ##########################################################
    time_Ma = timesteps.max() - timesteps
    w = correct_for_uranium_decay(time_Ma)

    if verbose is True:
        logger.info('w mean, min, max %0.2e, %0.2e, %0.2e' % (w.mean(), w.min(), w.max()))

    ###################################################
    # calculate observation frequency:
    ###################################################
    rho = calculate_normalized_density(rc)

    if verbose is True:
        logger.info('calculated observation frequency')
        logger.info('rho mean, min, max %0.2e, %0.2e, %0.2e' % (rho.mean(), rho.min(), rho.max()))

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
        logger.info('initial track length %0.2f' % l0)

    #################################
    # calculate track lengths
    #################################
    l = rm * l0
    # remove short track lengths (<2.18 um)
    l = np.where(l < min_length, 0, 1) * l

    # calculate track length standard deviation
    if apply_c_axis_correction is True:
        l_std = 0.008452 * l ** 2 - 0.2442 * l + 2.312
        l_std[np.where(l_std > 3.)] = 3.
    else:
        l_std = 0.02858 * l ** 2 - 0.8733 * l + 7.464
        l_std[np.where(l_std > 4.)] = 4.

    #######################################
    # calculate probability density function  (Item 5 — vectorized)
    #######################################
    if verbose is True:
        logger.info('start calculation of pdf of track lengths:')

    nbins = int(20 / binsize)
    bins_ = np.arange(0, 20, binsize)               # shape (nbins,)

    # Compute all per-timestep PDFs in one broadcast call:
    # l[:, np.newaxis]     shape (nsteps, 1)
    # bins_[np.newaxis, :] shape (1, nbins)
    # result               shape (nsteps, nbins)
    track_ln_prob = scipy.stats.norm.pdf(
        bins_[np.newaxis, :], l[:, np.newaxis], l_std[:, np.newaxis]
    )

    # Correct for uranium decay and observation frequency:
    track_ln_prob *= (w * rho)[:, np.newaxis]

    # Sum over timesteps and normalise:
    track_length_pdf = track_ln_prob.sum(axis=0)
    track_length_pdf = track_length_pdf / track_length_pdf.sum()

    if verbose is True:
        logger.info('done calculating probability density track lengths')

    ####################################################################
    # mean and standard deviation of model track length  (Item 5 — vectorized)
    ####################################################################
    if verbose is True:
        logger.info('calculate mean and std of track length from PDF')

    # bin centres: j*binsize + 0.5*binsize
    bin_centers = bins_ + 0.5 * binsize             # shape (nbins,)
    l_mean = np.dot(bin_centers, track_length_pdf)
    l_mean_std = np.sqrt(np.dot(track_length_pdf, (bins_ - l_mean) ** 2))
    l_median = track_length_pdf.argmax() * binsize

    if verbose is True:
        logger.info('mean track length  =  %0.2f,  std =  %0.2f,  median: %0.3f' % (l_mean, l_mean_std, l_median))
    if np.isnan(l_mean) is True:
        logger.warning('track length calculation failed')

    ###################################################
    # calculate observation frequency:
    ###################################################
    rho_s = r_standard
    dt = dts

    #####################
    # calculate AFT ages:
    #####################
    rc_mid = rc.copy()
    rc_mid[0] = rc[0]
    rc_mid[1:] = (rc[1:] + rc[:-1]) * 0.5

    rho_age = calculate_normalized_density(rc_mid) * w

    if verbose is True:
        logger.info('dt mean, min, max %0.2e, %0.2e, %0.2e' % (dt.mean() / Myr, dt.min() / Myr, dt.max() / Myr))
        logger.info('rho_age mean, min, max %0.2e, %0.2e, %0.2e' % (rho_age.mean(), rho_age.min(), rho_age.max()))
        logger.info('dt sum = %0.2f' % (dt.sum() / Myr))

    aft_age_uncorrected = np.dot(dt, rho_age)
    aft_age_corrected = aft_age_uncorrected / rho_s
    aft_age_myr = aft_age_corrected / Myr

    if verbose is True:
        logger.info('AFT age = %0.2f My, avg rho = %0.3f, rho standard = %s' % (aft_age_myr, rho_age.mean(), rho_s))

    if aft_age_corrected == 0:
        track_length_pdf[:] = 0
        l_mean = 0
        l_mean_std = 0

    if all_timesteps is True:
        return (track_length_pdf, aft_age_myr, l_mean, l_mean_std, rm, rc,
                rho_age, dt, aft_ages_myr)

    return track_length_pdf, aft_age_myr, l_mean, l_mean_std, rm, rc, rho_age, dt
