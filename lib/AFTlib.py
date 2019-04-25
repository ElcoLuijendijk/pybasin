'''
python apatite fission track library

AFT library:
- calculate fission track ages
- fission track age uncertainty & pdf
- c-axis projection track lengths
- conversion of Dpar between various etching methods

Elco Luijendijk 
elco.luijendijk at gmail.com

'''

import math, pdb, itertools

import numpy as np
import matplotlib.pyplot as pl
import scipy.stats
from scipy.stats import norm


__author__ = "Elco Luijendijk"
__copyright__ = "Copyright 2013, Elco Luijendijk"
__license__ = "GPL"
__version__ = "0.1"
__email__ = "elco.luijendijk at gmail.com"
__status__ = "Development"


def calculate_AFT_age_fluence(Ns, Ni, neutron_fluence, neutron_fluence_std,
                              lambda_d=1.551e-10, I=7.252e-3,
                              sigma_f=584.25e-28, lambda_f=8.51e-17,
                              g=0.5):
    
    """
    Calculate apatite fission track age & error 
    
    Equation valid for data that were not calculated using the zeta calibration 
    method (i.e., direct determination of neutron fluence)
    
    see: 
    Galbraith (1984), Mathematical Geology 16(7), 653-669
    Galbraith and Laslett (1985), Nucl. Tracks 10, 361-363
    
    Parameters
    ----------
    Ns : array
        number of spontaneous tracks
    Ni
        number of induced tracks
    neutron_fluence
        Neutron fluence (...)
    neutron_fluence_std
        standard devation of neutron fluence
    lambda_d
    
    I
    
    sigma_f
    
    lambda_f
    
    g
        geometry factor
    
    Returns
    -------
    FT_age : array
        calculated fission track age
    FT_age_min2SE : array
        calculated fission track age - 2*standard error of age
    FT_age_plus2SE : array
        calculated fission track age + 2*standard error of age
    
    """
    
    # sum the constants of the age equation 
    a = 1.0 / lambda_d
    b = lambda_d / lambda_f * I * sigma_f * g
    
    # calculate gamma:
    # eq. 2, Galbraith (1984)
    gamma = np.log( neutron_fluence * Ns / Ni )
    
    # calculate fission track age
    # eq. 3, Galbraith (1984)
    AFT_age = a * np.log(1.0 + b * np.exp(gamma)) / 1.0e6
    
    # calculate standard error of gamma
    # eq. 2 by Galbraith & Laslett (1985) Nucl. Tr. 10
    re_neutron_fluence = neutron_fluence_std / neutron_fluence
    se_gamma = np.sqrt(re_neutron_fluence **2 + 1.0/Ns + 1.0/Ni)

    # calculate FT age error:
    AFT_age_min2SE = a * np.log(1 + b * np.exp(gamma-2*se_gamma)) /1.0e6
    AFT_age_plus2SE = a * np.log(1 + b * np.exp(gamma+2*se_gamma)) / 1.0e6
    
    return AFT_age, AFT_age_min2SE, AFT_age_plus2SE


def calculate_AFT_age(Ns, Ni, Nd, rho_d, zeta, zeta_std):

    """
    Calculate apatite fission track ages using zeta calibration method
    
    equation 32 by Ketcham (2005)
    
    Parameters
    ----------
    Ns : array
    
    Ni : array
    
    Nd : array
    
    rho_d : array
    
    zeta : array
    
    zeta_std : array
    
    Returns
    -------
    AFT_age : array
    
    AFT_age_min1SE : array
    
    AFT_age_plus1SE : array
    
    AFT_age_min2SE : array
    
    AFT_age_plus2SE : array
    
    gamma : array
    
    gamma_std : array
    
    gamma_plus2SE : array
    """

    # calculate gamma:
    gamma = calculate_AFT_age_gamma(zeta, rho_d, Ns, Ni)
    gamma_std = calculate_AFT_age_gamma_std( gamma, Ns, Ni, Nd, 
                                                zeta, zeta_std )
    gamma_min2SE = gamma-2*gamma_std
    gamma_plus2SE = gamma+2*gamma_std
    gamma_min1SE = gamma-1*gamma_std
    gamma_plus1SE = gamma+1*gamma_std
    
    # calculate standard error of zero ages:
    Ns_dummy = np.ones(Ns.shape) * 3.0
    gamma_plus2SE[Ns==0] = calculate_AFT_age_gamma(zeta[Ns==0],
                                rho_d[Ns==0], Ns_dummy, Ni[Ns==0])
    Ns_dummy = np.ones(Ns.shape) * 1.5
    gamma_plus1SE[Ns==0] = calculate_AFT_age_gamma(zeta[Ns==0],
                                rho_d[Ns==0], Ns_dummy, Ni[Ns==0])
    
    # calculate AFT ages:
    AFT_age = calculate_AFT_age_from_gamma(gamma)
    AFT_age_min1SE = calculate_AFT_age_from_gamma(gamma_min1SE)
    AFT_age_plus1SE = calculate_AFT_age_from_gamma(gamma_plus1SE)
    AFT_age_min2SE = calculate_AFT_age_from_gamma(gamma_min2SE)
    AFT_age_plus2SE = calculate_AFT_age_from_gamma(gamma_plus2SE)
    
    # set zero ages for samples with 0 track counts:
    AFT_age[Ns==0] = 0.0
    AFT_age_min1SE[Ns==0] = 0.0
    AFT_age_min1SE[Ns==0] = 0.0
    
    return AFT_age, AFT_age_min1SE, AFT_age_plus1SE,\
            AFT_age_min2SE, AFT_age_plus2SE,\
            gamma, gamma_std, gamma_plus2SE


def calculate_AFT_age_gamma(zeta, rho_d, Ns, Ni):
    
    """
    Parameters
    ----------
    zeta : array or float
    
    rho_d : array or float
    
    Ns : array or float
        spontaneous tracks
    Ni : array or float
        induced tracks
        
    Returns
    -------
    gamma
    """
    
    if type(Ni) == np.ndarray:
        a = Ns/Ni.astype(float)
    else:
        a = Ns / float(Ni)
    gamma = np.log(zeta*rho_d*(a))
     
    return gamma


def calculate_AFT_age_SE(AFT_age, Ns, Ni, Nd, zeta, zeta_std):

    """
    Calculate standard error of fission track ages
    
    
    Parameters
    ----------
    AFT_age : array
        fission track age (Ma)
    Ns : array
        number of spontaneous tracks
    Ni
        number of induced tracks
    Nd
        
    zeta
        zeta
    zeta_std
        standard deviation of zeta
    
    Returns
    -------
    AFT_age_SE : array
        standard error of fission track age
    """
    
    # calculate standard error of zeta:
    SE_zeta = zeta_std/zeta
    
    # calculate standard error of gamma:
    AFT_age_SE = (AFT_age * 1.0e6 * np.sqrt((Ns**-1) + (Ni**-1) +
                    (Nd**-1) + (SE_zeta**2))) / 1.0e6
    
    # calculate standard error gamma for samples with 0 counted tracks:
    ind_ = np.where(Ns == 0)[0]
    if len(ind_)>0:
        Ns_dummy = np.ones(len(ind_)) * 3.
        AFT_age_SE[ind_] = (AFT_age[ind_] * 1.0e6 *
                            np.sqrt((Ns_dummy**-1) + (Ni[ind_]**-1) +
                            (Nd[ind_]**-1)+(SE_zeta[ind_]**2))) / 1.0e6
        
    return AFT_age_SE


def calculate_AFT_age_gamma_std(gamma, Ns, Ni, Nd, zeta, zeta_std):
    
    """
    Parameters
    ----------
    gamma : array
    
    Ns : array
        number of spontaneous tracks
    Ni : array
        number of induced tracks    
    Nd : array
    
    zeta : array
    
    zeta_std : array
    
    
    Returns
    -------
    gamma_std : array

    """
    
    # calculate standard error of zeta:
    SE_zeta = zeta_std/zeta
    
    # calculate standard error of gamma:
    gamma_std = np.sqrt((Ns**-1)+(Ni**-1)+(Nd**-1)+(SE_zeta)**2)
    
    ## calculate standard error gamma for samples with 0 counted tracks:
    # check if input is numpy array:
    if type(Ns) == np.ndarray:
        Ns_dummy = 3
        gamma_std_zero = np.sqrt(   Ns_dummy**-1 + Ni**-1 + 
                                    Nd**-1 + (SE_zeta**2) )
        gamma_std = np.where(Ns == 0, gamma_std_zero, gamma_std)

    else:
        # non-array input:
        Ns_dummy = 3
        gamma_std = np.sqrt((Ns_dummy**-1) +
                                            (Ni**-1) +
                                            (Nd**-1) +
                                            (SE_zeta**2))

    return gamma_std


def calculate_AFT_age_from_gamma(gamma):
    
    """
    FT age equation
    
    recast following    Galbraith (1984) Math. Geol. 16(7) and 
                        Galbraith&Laslett (1985) Nuclear Tracks 10
    
    Parameters
    ----------
    gamma : array
    
    Returns
    -------
    age : array
    
    """
    
    
    lambda_d = 1.551e-10
    a =  1./lambda_d
    b = 0.5*lambda_d   # geometry factor 
    
    age = a*np.log(1+b*np.exp(gamma))/1.0e6
    
    return age

    
def calculate_gamma_from_AFT_age(age, lambda_d = 1.551e-10, g=0.5):
    
    """
    Parameters
    ----------
    age : array
    
    lambda_d = 1.551e-10 : float
    
    g=0.5 : float
    
    Returns
    -------  
    gamma  : array
    """
    
    
    a =  1./lambda_d
    b = g*lambda_d   # geometry factor 
    
    #age = a*np.log(1+b*np.exp(gamma))/1.0e6
    gamma = np.log((np.exp(age*1.0e6/a)-1) / b)
    
    return gamma
  

def calculate_central_age( Ns, Ni, Nd, rho_d, zeta, zeta_se,
                           lambda_d=1.551e-10, g= 0.5,
                           returnEta = False):
    
    '''
    Algorithm for calculating central fission track ages
    
    See: Galbraith & Laslett (1993) International Journal of Radiation
    Applications and Instrumentation. 21, 4: 459-470.
    doi:10.1016/1359-0189(93)90185-C.
  
    Parameters
    ----------
    Ns : array
    
    Ni : array
    
    Nd : array
    
    rho_d : array
    
    zeta : array
    
    zeta_se : array
    
    lambda_d=1.551e-10 : float
    
    g= 0.5 : float
        geometry factor
    
    Returns
    -------    
    central_age : array
        central fission track age
    central_age_se : array
        standard error of central fission track age
    
    '''
    
    Nsamples = len(Ns)
    
    Ns = Ns.astype(float)
    Ni = Ni.astype(float)
    Nd = Nd.astype(float)
    
    m = Ns+Ni
    y = Ns/m
    z = np.log(Ns+0.5) - np.log(Ni+0.5)
    
    sigma = 0.6 * z.std()
    eta = Ns.sum()/m.sum()
    
    u = np.zeros((Nsamples))
    w = np.zeros((Nsamples))
    
    Niterations = 0
    diff_sigma = 1000.0
    diff_eta = 1000.0
    
    while Niterations < 20 and abs(diff_sigma) >1.0e-3 and abs(diff_eta)>1.0e-3:
        for i, mi, yi in zip(itertools.count(), m, y):
            
            w[i] = mi / (eta*(1-eta)+(mi-1)*eta**2*(1-eta)**2*sigma**2)
            
            u[i] = w[i]**2*(yi-eta)**2
        
        diff_sigma = ((u.sum()/w.sum())**0.5) - 1.0
        sigma = sigma * ((u.sum()/w.sum())**0.5)
        
        new_eta = (w*y).sum() / w.sum()
        diff_eta = eta/new_eta  -1.0
        
        eta = new_eta

        Niterations += 1
        
        if Niterations>100:
            print('error, too many iterations')
            pdb.set_trace()
    
    central_age = 1.0/lambda_d * np.log(1.0 + lambda_d * g * zeta * rho_d.mean() *(eta/(1.0-eta))) /1.0e6
    
    rse = (1.0 / (eta **2 * (1-eta)**2 *w.sum()) + (1.0/Nd.sum()) + (zeta_se/zeta)**2) ** 0.5
    
    central_age_se = rse * central_age

    if returnEta == True:
        return eta
    else:
        return central_age, central_age_se


def chi_sq_test(Nsj, Nij):
    
    '''
    Chi square test of fission track samples
    
    see eq. 8 and 9 of Green (1981), Nuclear Tracks 5(1-2) p. 77-86

    
    Parameters
    ----------
    Nsj : array
        number of spontaneous tracks
    Nij : array or float
        number of induced tracks
        
    Returns
    -------
    p_chisq : float
        chi squared probability
        
    '''
    
    
    Ns = float(Nsj.sum())
    Ni = float(Nij.sum())
    Ngrains = len(Nsj)
    chi_sq = 0
    
    for j in range(Ngrains):
        
        Nsj_e = Ns / (Ns + Ni) * (Nsj[j] + Nij[j])
        Nij_e = Ni / (Ns + Ni) * (Nsj[j] + Nij[j])
        
        chi_sq += (Nsj[j] - Nsj_e) **2 / Nsj_e + (Nij[j] - Nij_e) **2 / Nij_e
    
    dof = Ngrains - 1

    p_chisq = 1.0 - scipy.stats.chi2.cdf(chi_sq, dof)
    
    return p_chisq


def APFU_to_Cl_wt_fraction(Cl_apfu): #, F_apfu, OH_apfu):
    
    """
    convert Chloride APFU units to weight fraction

    Parameters
    ----------
    Cl_apfu : array or float

    Returns
    -------  
    Cl_wtfract : array or float

    """
    
    # atomic weights (g mol-1)
    Ca = 40.078
    P = 30.973762
    O = 15.9994
    H = 1.00794
    F = 18.9984032
    Cl = 35.453
    Br = 79.904

    # weight fluorapatite:
    #apatite_weight = Ca*10 + (P+4*O)*6 + F*F_apfu + Cl*Cl_apfu + (O+H)*OH_apfu
    apatite_weight = Ca*10 + (P+4*O)*6 + F*2 #+ Cl*Cl_apfu + (O+H)*OH_apfu
    Cl_weight = Cl*Cl_apfu
    Cl_wtfract = Cl_weight / apatite_weight

    return Cl_wtfract
    

def Cl_wt_fraction_to_APFU(Cl_wtfract):
    
    '''
    convert Chloride weight fraction to APFU units
    
    Parameters
    ----------
    Cl_wtfract
        weight fraction of chloride in apatite
    
    Returns
    -------  
    Cl_apfu
        chloride in apfu units

    '''
    
    # atomic weights (g mol-1)
    Ca = 40.078
    P = 30.973762
    O = 15.9994
    H = 1.00794
    F = 18.9984032
    Cl = 35.453
    Br = 79.904

    # weight fluorapatite:
    apatite_weight_noAnions = Ca*5 + (P+4*O)*3
    
    Cl_apfu = (-Cl_wtfract*apatite_weight_noAnions - Cl_wtfract*F) / \
        (Cl_wtfract*Cl - Cl - Cl_wtfract*F)
    

    return Cl_apfu


def convert_Dpar_Barbarand2003(Dpar50):
    
    '''
    
    Convert Dpar from  Barbarand et al (2003) 5.0M etching conditions
    to Carlson et al (1999) 5.5 M etching conditions

    See Ketcham et al. (2007) Am. Min. 92: 799-810
    
    Parameters
    ----------
    Dpar50 : array
        Dpar measured following Barbarand et al.(2003) 5.0M HNO3, .. sec etching
        coditions (um)
    
    Returns
    -------  
    Dpar55 : array
        Dpar measured following Carlson et al.(1999) 5.5M HNO3, .. sec etching
        coditions (um)
    
    '''
    
    Dpar55 = 0.9231*Dpar50 + 0.2515
    
    return Dpar55


def convert_Dpar_VU(Dpar_VU, slope=1.12, slope_std=0.21 ):
    
    '''
    convert Dpar values of the Vrije Universiteit (VU) Amsterdam
    lab to equivalent Donelick (1999) values.
    
    etching methods
    - VU: 1.6M (7%) HNO3, 35 sec
    - Donelick: 5.5M HNO3, 20 sec
    
    based on linear correlation using data from:
    Murrell (2009). Geol. Soc. Spec. Pub. 1, 73-85
    
    takes into account +-95% confidence interval of correlation
    
    Parameters
    ----------
    Dpar_VU : array or float
        Dpar (um) measured using VU etching methods, 1.6M HNO3, 35 sec 
    slope=1.12 : float
        slope regression line between VU and Donelick etching methods
    slope_std=0.21 : float
        standard devation of slope regression line VU and Donelick methods
        
    Returns
    -------  
    Dpar_Donelick_min : array or float
        min estimate of Dpar (um)
    Dpar_Donelick_max : array or float
        max estimate of Dpar (um)
    
    '''

    Dpar_Donelick_min = (slope - slope_std) * Dpar_VU_min 
    Dpar_Donelick_max = (slope + slope_std) * Dpar_VU_max

    return Dpar_Donelick_min, Dpar_Donelick_max
    

def caxis_proj_reduced_lengths(r, a=-1.499, b=4.150, c=-1.656):
    
    '''
    Convert reduced lengths to c-axis projected reduced lengths
    
    See Fig 8 in Ketcham et al. (1999) Am. Mineralogist 84, 1235-1255.
    
    Parameters
    ----------
    r : array or float
        reduced track length (um)
    a= -1.499 : float
        constant
    b=4.150 : float
        constant
    c=-1.656 : float
        constant
        
    Returns
    -------  
    rm : array or float
        c-axis projected reduced track length (um)
    
    '''
    
    rm = a*(r**2) + b*r + c
    
    return rm


def caxis_proj_lengths(l, a=-0.08076, b = 3.856, c=-25.488):
    
    '''
    Convert track lengths to c-axis projected lengths
    
    see Fig 8 in Ketcham et al. (1999) Am. Mineralogist 84, 1235-1255.
    

    Parameters
    ----------
    l : array or float
        track length (um)
    a=-0.08076 : float
        constant
    b = 3.856 : float
        constant
    c=-25.488 : float
        constant
    
    Returns
    -------  
    lm : array_like
        c-axis projected track length

    '''
    
    cl = c-l
    D = (b**2-(4*a*cl))
    if D<0:
        print('error,  D<0')
        print(bla)
        
    lm = (-b+(math.sqrt(D)))/(2*a)
    if lm<0:
        lm = (-b-(math.sqrt(D)))/(2*a)
    
    return lm

   
def calculate_AFT_ages_PDF(bins, AFT_ages, AFT_ages_min, AFT_ages_max):
    
    '''
    calculate probability density function of AFT age data
    
    based on:
    - Galbraith, R F. 1984. Mathematical Geology 16, no. 7: 653-669.
    - Galbraith, R F, G. M. Laslett. 1985. Nuclear Tracks 10: 361-363.
    
    Parameters
    ----------
    bins : numpy array
        bins of the probability density function
    AFT_ages : numpy array
        fission track ages (Ma?)
    AFT_ages_min : numpy array
        ?
    AFT_ages_max : numpy array
        ?
    
    Returns
    -------  
    gamma_pdf : numpy array
        probability density function of gamma
    
    '''
    
    AFT_ages_min_ = AFT_ages+AFT_ages_min
    AFT_ages_max_ = AFT_ages+AFT_ages_max
    
    # calculate gamma value (i.e. normally distributed transform of AFT age equation)
    gamma = calculate_gamma_from_AFT_age(AFT_ages)
    
    # correct zero ages:
    gamma[AFT_ages == 0] = 1e-10
    
    # calculate standard error of gamma:
    gamma_max = calculate_gamma_from_AFT_age(AFT_ages_max_)
    gamma_std = (gamma_max - gamma) / 2.0
    
    # calculate pdf of gamma:
    gamma_pdf = np.zeros(bins.shape)
    gamma_bins = calculate_gamma_from_AFT_age(bins)
    for gamma_i, gamma_std_i in zip( gamma, gamma_std ) :
        gamma_pdf += norm.pdf(gamma_bins, gamma_i, gamma_std_i)
    
    # normalize pdf to a value of 1:
    gamma_pdf = gamma_pdf / gamma_pdf.sum()
    
    return gamma_pdf
