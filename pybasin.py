"""
PyBasin, version 0.1

recoded in 2014-2015 to simplify the code and to use .csv files
+ pandas for model input

Elco Luijendijk, Goettingen University

<elco.luijendijk@geo.uni-goettingen.de>

"""

import sys
import os
import argparse
import imp
import ast
# from runpy import run_path


# import pdb
import datetime
import pickle
import itertools
import time
import inspect
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as pl

from multiprocessing import Pool

import lib.pybasin_lib as pybasin_lib
import lib.pybasin_figures as pybasin_figures

# helium diffusion algortihm by Meesters and Dunai (2003)
try:
    import lib.helium_diffusion_models as he
except ImportError:
    print('warning, failed to import AHe module')


# make sure multi-threading for numpy is turned off (this slows down the heat
# flow solution a lot...)
os.environ['OPENBLAS_NUM_THREADS'] = '1'


def model_data_comparison_T(T_data_well, z_nodes, T_nodes, active_nodes):

    """
    Compare modelled and measured subsurface temperatures

    """

    T_data_well['simulated_T'] = np.interp(T_data_well['depth'],
                                           z_nodes[-1, active_nodes[-1]],
                                           T_nodes[-1, active_nodes[-1]])

    # calculate model error temperature data
    # ind = [T_data_well['temperature_unc_1sigma'].isnull()]
    # T_data_well['temperature_unc_1sigma'][ind] = pybasin_params.vr_unc_sigma
    T_data_well['residual'] = (T_data_well['temperature']
                               - T_data_well['simulated_T'])
    T_data_well['residual_norm'] = (T_data_well['residual']
                                    / T_data_well['temperature_unc_1sigma'])
    T_data_well['P_fit'] = \
        (1.0 - scipy.stats.norm.cdf(np.abs(T_data_well['residual_norm']))) * 2

    # assign P=0 for temperatures lower than uncorrected BHTs
    ind_bht_nofit = ((T_data_well['residual'] > 0)
                     & (T_data_well['data_type'] == 'BHT'))
    ind_bht_ok = ((T_data_well['residual'] <= 0)
                  & (T_data_well['data_type'] == 'BHT'))
    T_data_well['P_fit'][ind_bht_nofit] = 0
    T_data_well['residual'][ind_bht_nofit] = 15.0
    T_data_well['P_fit'][ind_bht_ok] = 1.00
    T_data_well['residual'][ind_bht_ok] = 0.0

    T_rmse = np.sqrt(np.mean(T_data_well['residual']**2))
    T_gof = np.mean(T_data_well['P_fit'])

    print('Temperature data:')
    print(T_data_well)

    return T_gof, T_rmse


def model_data_comparison_VR(vr_data_well, z_nodes, vr_nodes, active_nodes, vr_unc_sigma=0.05):

    """
    Compare modeled and measured vitrinite reflectance values
    """

    vr_data_well['simulated_vr'] = \
        np.interp(vr_data_well['depth'],
                  z_nodes[-1, active_nodes[-1]],
                  vr_nodes[-1, active_nodes[-1]])

    # calculate model error vitrinite data
    vr_data_well['residual'] = (vr_data_well['VR']
                                - vr_data_well['simulated_vr'])

    # if min and max VR is given calculate asymetric error value:
    indm = vr_data_well['VR_min'].notnull() & vr_data_well['VR_max'].notnull()
    vr_data_well['vr_SE_plus'] = vr_data_well['VR_max'] - vr_data_well['VR']
    vr_data_well['vr_SE_min'] = vr_data_well['VR'] - vr_data_well['VR_min']
    ind_plus = indm & vr_data_well['residual'] >= 0
    ind_neg = indm & vr_data_well['residual'] < 0

    # min and max value assumed to be +- 2 SE
    vr_data_well.loc[ind_plus, 'VR_unc_1sigma'] = \
        vr_data_well.loc[ind_plus, 'vr_SE_plus'] / 2.0
    vr_data_well.loc[ind_neg, 'VR_unc_1sigma'] = \
        vr_data_well.loc[ind_neg, 'vr_SE_min'] / 2.0

    # otherwise, use value of 1 sigma unc. given in input file:
    ind = (vr_data_well['VR_unc_1sigma'].isnull()) & (indm == False)
    vr_data_well.loc[ind, 'VR_unc_1sigma'] = vr_unc_sigma

    # calculate normalized residual
    vr_data_well['residual_norm'] = (vr_data_well['residual'] / vr_data_well['VR_unc_1sigma'])
    vr_data_well['P_fit'] = (1.0 - scipy.stats.norm.cdf(np.abs(vr_data_well['residual_norm'])))*2

    # calculate total rmse and goodness of fit statistic:
    vr_rmse = np.sqrt(np.mean(vr_data_well['residual']**2))
    vr_gof = np.mean(vr_data_well['P_fit'])

    return vr_rmse, vr_gof, vr_data_well


def model_data_comparison_AFT_age(aft_data_well, aft_ages,
                                  modeled_aft_age_samples_min,
                                  modeled_aft_age_samples_max,
                                  verbose=False):

    """
    Compare modelled and measured apatite fission track ages
    """

    age_bins = []
    age_pdfs = []
    single_grain_aft_ages = []
    single_grain_aft_ages_se_min = []
    single_grain_aft_ages_se_plus = []

    # go through all samples
    for sample in aft_data_well['sample']:

        if sample in aft_ages['sample'].values:
            ind_sample = aft_ages['sample'].values == sample

            # find single grain ages for this sample
            single_grain_ages_sample = aft_ages['aft_age'][ind_sample].values
            single_grain_ages_se_min_sample = \
                aft_ages['aft_age_stderr_min'][ind_sample].values
            single_grain_ages_se_plus_sample = \
                aft_ages['aft_age_stderr_plus'][ind_sample].values
            single_grain_aft_ages.append(single_grain_ages_sample)
            single_grain_aft_ages_se_min.append(
                single_grain_ages_se_min_sample)
            single_grain_aft_ages_se_plus.append(
                single_grain_ages_se_plus_sample)

            # get pdf of observed AFT ages from single grain ages
            age_bin, age_pdf = \
                pybasin_lib.calculate_aft_ages_pdf(
                    single_grain_ages_sample,
                    single_grain_ages_se_min_sample,
                    single_grain_ages_se_plus_sample)

        else:
            # get pdf of observed age from central age instead
            ind_sample = aft_data_well['sample'].values == sample

            # get pdf of observed AFT ages
            age_bin, age_pdf = \
                pybasin_lib.calculate_aft_ages_pdf(
                    aft_data_well['aft_age'][ind_sample].values,
                    aft_data_well['aft_age_stderr_min'][ind_sample].values,
                    aft_data_well['aft_age_stderr_plus'][ind_sample].values)

            single_grain_aft_ages.append(None)
            single_grain_aft_ages_se_min.append(None)
            single_grain_aft_ages_se_plus.append(None)

            # calculate error for completely reset samples
            # (ie 0 age with SE of 0)
            if aft_data_well.loc[ind_sample, 'aft_age'].values == 0 \
                    and aft_data_well.loc[ind_sample, 'aft_age_stderr_plus'].values == 0:

                age_bin = np.array([0, 1e-3, 1e3])
                age_pdf = np.array([1.0, 0.0, 0.0])

        age_bins.append(age_bin)
        age_pdfs.append(age_pdf)

    # go through samples and find out how much of age pdf is covered by
    #  min and max simulated age
    for i, sample_ix, age_bin, age_pdf in zip(itertools.count(),
                                              aft_data_well.index,
                                              age_bins,
                                              age_pdfs):

        if np.any(np.isnan(age_pdf)) == False:
            #
            aft_data_well.loc[sample_ix, 'simulated_AFT_min'] = modeled_aft_age_samples_min[i]
            aft_data_well.loc[sample_ix, 'simulated_AFT_max'] = modeled_aft_age_samples_max[i]

            # TODO: find more elegant solution for 0.0 simulated AFT age
            # and check if GOF for AFT ages of 0.0 Ma are correct
            if aft_data_well.loc[sample_ix, 'simulated_AFT_min'] == 0:
                start_ind = 0
            else:
                start_ind = np.where(
                    aft_data_well.loc[sample_ix, 'simulated_AFT_min']
                    >= age_bin)[0][-1]

            # if aft_data_well.loc[sample_ix, 'simulated_AFT_max'] == 0.0:
            #    end_ind = 0
            # else:
            #    np.where(0.0 >= age_bins)[0]
            end_ind = np.where(aft_data_well.loc[sample_ix, 'simulated_AFT_max'] < age_bin)[0][0]

            pdf_fit_sum = np.sum(age_pdf[start_ind:end_ind])
            pdf_nofit_sum = np.sum(age_pdf[:start_ind]) + np.sum(age_pdf[end_ind:])
            aft_data_well.loc[sample_ix, 'GOF_aft_ages'] = pdf_fit_sum

            # if aft_data_well.loc[sample_ix, 'aft_age'] == 0 \
            #        and aft_data_well.loc[sample_ix, 'aft_age_stderr_plus'] == 0:

            #    print modeled_aft_age_samples_min[i],  modeled_aft_age_samples_max[i], pdf_fit_sum
        else:
            aft_data_well.loc[sample_ix, 'GOF_aft_ages'] = np.nan

    # calculate model error:
    for i, sample_ix, age_bin, age_pdf in zip(itertools.count(),
                                              aft_data_well.index,
                                              age_bins,
                                              age_pdfs):

        if np.any(np.isnan(age_pdf)) == False:

            pc = np.cumsum(age_pdf)

            if pc[0] == 1:
                start_ind = 0
                end_ind = 1
            else:
                # find +-95% confines of age distribution
                start_ind = np.where(pc >= 0.05)[0][0]
                end_ind = np.where(pc <= 0.95)[0][-1]

            age_min = age_bin[start_ind]
            age_max = age_bin[end_ind]

            # check difference of min modeled aft age and min. value of age
            # distribution
            if modeled_aft_age_samples_min[i] < age_min:
                age_error_min = 0
            else:
                age_error_min = modeled_aft_age_samples_min[i] - age_min

            # check difference of max modeled aft age and max. value of age
            # distribution
            if modeled_aft_age_samples_max[i] > age_max:
                age_error_max = 0
            else:
                age_error_max = age_max - modeled_aft_age_samples_max[i]

            # differerent procedure for observed AFT ages of 0, add penaly
            # if modeled ages are older:
            # check difference of min modeled aft age and min. value of age
            # distribution
            if age_max <= 1e-3:
                age_error_max = modeled_aft_age_samples_max[i]

            age_error = age_error_min + age_error_max

            aft_data_well.loc[sample_ix, 'age_error'] = age_error

        else:
            aft_data_well.loc[sample_ix, 'age_error'] = np.nan
            print('no model-data comparison for sample %s, '
                  'missing age data?' % sample_ix)

    # calculate mean GOF from single grain GOFs for each sample
    aft_age_mean_gof = aft_data_well['GOF_aft_ages'].dropna().mean()
    aft_age_mean_error = aft_data_well['age_error'].dropna().mean()

    if verbose is True:

        print(aft_data_well[['sample', 'depth',
                            'aft_age', 'aft_age_stderr_plus',
                            'simulated_AFT_min', 'simulated_AFT_max',
                            'GOF_aft_ages', 'age_error']])

    return (aft_age_mean_gof, aft_age_mean_error,
            single_grain_aft_ages, single_grain_aft_ages_se_min,
            single_grain_aft_ages_se_plus,
            age_bins, age_pdfs, aft_data_well)


def model_data_comparison_AHe(ahe_samples_well, ahe_data,
                              ahe_age_bin,
                              modeled_ahe_age_samples_min,
                              modeled_ahe_age_samples_max):

    """
    Compare modelled and measure apatite (U-Th)/He ages

    """

    print('calculating GOF AHe data')

    ahe_age_pdfs_all_samples = []
    ahe_ages_all_samples = []
    ahe_ages_all_samples_SE = []

    for ahe_sample_i, ahe_sample_ix, ahe_sample in zip(
            itertools.count(),
            ahe_samples_well.index,
            ahe_samples_well['sample']):

        ahe_age_pdfs = []

        if ahe_sample in ahe_data['sample'].values:
            ind_sample = ahe_data['sample'].values == ahe_sample

            grain_pdfs = []

            ahe_ages_all_samples.append(
                ahe_data.loc[ind_sample, 'ahe_age_uncorrected'].values)
            ahe_ages_all_samples_SE.append(
                ahe_data.loc[ind_sample, 'ahe_age_uncorrected_se'].values)

            age_error = 0

            for grain_i, ahe_age_obs, ahe_age_obs_SE \
                    in zip(itertools.count(),
                           ahe_data.loc[ind_sample, 'ahe_age_uncorrected'].values,
                           ahe_data.loc[ind_sample, 'ahe_age_uncorrected_se'].values):

                ahe_age_pdf = scipy.stats.norm.pdf(ahe_age_bin,
                                                   ahe_age_obs,
                                                   ahe_age_obs_SE)

                # normalize to make sum of pdf 1
                ahe_age_pdf = ahe_age_pdf / ahe_age_pdf.sum()

                ahe_age_pdfs.append(ahe_age_pdf)

                # find out how much of pdf is covered by simulated
                # end-member AHe ages
                # ahe_age_pdf = ahe_age_pdf / ahe_age_pdf.sum()

                ahe_age_sim_min = \
                    modeled_ahe_age_samples_min[ahe_sample_i][grain_i]
                ahe_age_sim_max = \
                    modeled_ahe_age_samples_max[ahe_sample_i][grain_i]

                if ahe_age_sim_min == 0:
                    start_ind = 0
                else:
                    start_ind = np.where(ahe_age_sim_min >= ahe_age_bin)[0][-1]

                if ahe_age_sim_max == 0.0:
                    end_ind = 0
                else:
                    end_ind = np.where(ahe_age_sim_max <= ahe_age_bin)[0][0]

                pdf_fit_sum = np.sum(ahe_age_pdf[start_ind:end_ind])

                grain_pdfs.append(pdf_fit_sum)

                # calculate ahe age error:
                # if np.any(np.isnan(age_pdf)) == False:
                # pc = np.cumsum(ahe_age_pdf)

                # find +-95% confines of age distribution
                # start_ind = np.where(pc >= 0.05)[0][0]
                # end_ind = np.where(pc <= 0.95)[0][-1]

                age_min = ahe_age_obs - ahe_age_obs_SE * 1.96
                age_max = ahe_age_obs + ahe_age_obs_SE * 1.96

                # check difference of min modeled aft age and min. value of age distribution
                if ahe_age_sim_min < age_min:
                    age_error_min = 0
                else:
                    age_error_min = ahe_age_sim_min - age_min

                # check difference of max modeled aft age and max. value of age distribution
                if ahe_age_sim_max > age_max:
                    age_error_max = 0
                else:
                    age_error_max = age_max - ahe_age_sim_max

                age_error += age_error_min + age_error_max
                # aft_data_well.loc[sample_ix, 'age_error'] = age_error

            ahe_samples_well.loc[ahe_sample_ix, 'mean_GOF_all_grains'] = \
                np.mean(np.array(grain_pdfs))
            ahe_samples_well.loc[ahe_sample_ix, 'min_GOF_all_grains'] = \
                np.min(np.array(grain_pdfs))
            ahe_samples_well.loc[ahe_sample_ix, 'max_GOF_all_grains'] = \
                np.max(np.array(grain_pdfs))
            ahe_samples_well.loc[ahe_sample_ix, 'mean_ahe_error'] = age_error / len(grain_pdfs)

        ahe_age_pdfs_all_samples.append(ahe_age_pdfs)

    if 'mean_GOF_all_grains' in ahe_samples_well.columns:
        ahe_age_gof = ahe_samples_well['mean_GOF_all_grains'].mean()
    else:
        ahe_age_gof = np.nan

    if 'mean_ahe_error' in ahe_samples_well.columns:
        ahe_age_error = ahe_samples_well.loc[ahe_sample_ix, 'mean_ahe_error'].mean()
    else:
        ahe_age_error = 99999.9

    return (ahe_age_gof, ahe_age_error, ahe_ages_all_samples, ahe_ages_all_samples_SE,
            ahe_age_bin, ahe_age_pdfs_all_samples, ahe_samples_well)


def model_data_comparison_salinity(salinity_data_well,
                                   z_nodes, C_nodes, active_nodes):

    """
    Compare modeled and measured porewater salinity
    """

    salinity_data_well['simulated_salinity'] = \
        np.interp(salinity_data_well['depth'],
                  z_nodes[-1, active_nodes[-1]],
                  C_nodes[-1, active_nodes[-1]])

    # calculate model error salinity data
    salinity_data_well['residual'] = \
        (salinity_data_well['salinity']
         - salinity_data_well['simulated_salinity'])
    salinity_data_well['residual_norm'] = \
        (salinity_data_well['residual']
         / salinity_data_well['salinity_unc_1sigma'])
    salinity_data_well['P_fit'] = \
        (1.0 - scipy.stats.norm.cdf(
            np.abs(salinity_data_well['residual_norm']))) * 2

    salinity_rmse = \
        np.sqrt(np.mean(salinity_data_well['residual']**2))
    salinity_gof = np.mean(salinity_data_well['P_fit'])

    return salinity_gof


def assemble_data_and_simulate_aft(resample_t, nt_prov,
                                   n_nodes, time_array_bp,
                                   z_nodes, T_nodes, active_nodes,
                                   prov_start_nodes, prov_end_nodes,
                                   annealing_kinetics_values,
                                   annealing_kinetic_param,
                                   surface_temp,
                                   aft_data_well,
                                   calculate_thermochron_for_all_nodes=False,
                                   annealing_eq='FC',
                                   C0=0.39528, C1=0.01073,
                                   C2=-65.12969, C3=-7.91715,
                                   alpha=0.04672,
                                   location_has_AFT=True):

    """
    Use modeled temperature history and provneance history to model apatite fission track ages and length distributions
    """

    # resample_t = pybasin_params.resample_timesteps
    # nt_prov = pybasin_params.provenance_time_nt
    # pybasin_params.annealing_kinetics_values,
    # pybasin_params.annealing_kinetic_param,

    # pybasin_params.make_model_data_fig is True:

    if calculate_thermochron_for_all_nodes is True:
        print('calculating AFT ages and lengths for all n=%i nodes' % n_nodes)
        # simulate AFT all nodes
        simulated_AFT_data =\
            pybasin_lib.simulate_aft(
                resample_t, nt_prov, n_nodes, time_array_bp,
                z_nodes, T_nodes, active_nodes,
                prov_start_nodes, prov_end_nodes,
                annealing_kinetics_values,
                annealing_kinetic_param,
                surface_temp,
                annealing_eq=annealing_eq,
                C0=C0, C1=C1, C2=C2, C3=C3, alpha=alpha)

        (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
         aft_ln_mean_nodes, aft_ln_std_nodes,
         aft_node_times_burial, aft_node_zs,
         aft_node_times, aft_node_temps) = \
            simulated_AFT_data

    else:
        simulated_AFT_data = None

    nt = T_nodes.shape[0]
    n_aft_samples = len(aft_data_well)

    if n_aft_samples == 0:
        return (None, None, None, None, None, None, None, None, None,
                simulated_AFT_data, None, None)

    # get T history for samples only
    T_samples = np.zeros((nt, n_aft_samples))
    for h in range(nt):
        T_samples[h, :] = \
            np.interp(aft_data_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      T_nodes[h, active_nodes[-1]])

    # get burial history of samples
    z_aft_samples = np.zeros((nt, n_aft_samples))
    for h in range(nt):
        z_aft_samples[h, :] = \
            np.interp(aft_data_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      z_nodes[h, active_nodes[-1]])

    # get provenance history for samples only
    n_prov = prov_start_nodes.shape[1]
    prov_start_samples = np.zeros((n_aft_samples, n_prov))
    prov_end_samples = np.zeros((n_aft_samples, n_prov))
    for h in range(n_prov):
        prov_start_samples[:, h] = \
            np.interp(aft_data_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      prov_start_nodes[active_nodes[-1], h])
        prov_end_samples[:, h] = \
            np.interp(aft_data_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      prov_end_nodes[active_nodes[-1], h])

    # get active node array for samples only
    active_nodes_aft_samples = np.zeros((nt, n_aft_samples),
                                        dtype=bool)
    for h in range(nt):
        active_nodes_aft_samples[h, :] = \
            np.interp(aft_data_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      active_nodes[h, active_nodes[-1]])

    # select prov. history for samples:
    print('calculating AFT ages and lengths for n=%i samples' % n_aft_samples)
    simulated_aft_data_samples =\
        pybasin_lib.simulate_aft(
            resample_t, nt_prov, n_aft_samples, time_array_bp,
            z_aft_samples, T_samples, active_nodes_aft_samples,
            prov_start_samples, prov_end_samples,
            annealing_kinetics_values,
            annealing_kinetic_param,
            surface_temp,
            annealing_eq=annealing_eq,
            C0=C0, C1=C1, C2=C2, C3=C3, alpha=alpha)

    (modeled_aft_age_samples, modeled_aft_age_samples_min,
     modeled_aft_age_samples_max,
     aft_ln_mean_samples, aft_ln_std_samples,
     aft_sample_times_burial, aft_sample_zs,
     aft_sample_times, aft_sample_temps) = simulated_aft_data_samples

    return (modeled_aft_age_samples,
            modeled_aft_age_samples_min,
            modeled_aft_age_samples_max,
            aft_ln_mean_samples,
            aft_ln_std_samples,
            aft_sample_times_burial,
            aft_sample_zs,
            aft_sample_times,
            aft_sample_temps,
            simulated_AFT_data,
            z_aft_samples,
            T_samples)


def assemble_data_and_simulate_AHe(ahe_samples_well,
                                   ahe_data,
                                   decay_constant_238U,
                                   decay_constant_235U,
                                   decay_constant_232Th,
                                   n_nodes, resample_t, nt_prov,
                                   time_array_bp,
                                   z_nodes,
                                   T_nodes,
                                   active_nodes,
                                   prov_start_nodes,
                                   prov_end_nodes,
                                   surface_temp,
                                   calculate_thermochron_for_all_nodes=False,
                                   U_default=50.0, Th_default=50.0,
                                   radius_default=75.0,
                                   ahe_method='RDAAM',
                                   alpha=0.04672, C0=0.39528, C1=0.01073,
                                   C2=-65.12969, C3=-7.91715):

    """
    Use modeled temperature history and provneance history to model apatite (U-Th)/He ages
    """

    if calculate_thermochron_for_all_nodes is True:

        print('-' * 10)
        print('calculating AHe ages for all nodes')

        ahe_grain_radius_nodes = np.zeros((n_nodes, 2))
        U_nodes = np.zeros((n_nodes, 2))
        Th_nodes = np.zeros((n_nodes, 2))

        # fill with default values
        # TODO: specify default values in input file
        U_nodes[:, :] = U_default
        Th_nodes[:, :] = Th_default
        ahe_grain_radius_nodes[:, :] = radius_default

        Ur0_max = 0
        Ur0_min = 99999

        # find min and max grain diameters and U and Th contents
        # for this location
        samples = ahe_samples_well['sample'].values
        # ahe_grain_radius_samples = []
        for ahe_sample_no, ahe_sample in enumerate(samples):
            ind_sample = ahe_data['sample'] == ahe_sample
            ahe_grain_radius_sample = \
                ahe_data['grain_radius'][ind_sample].values * 1e-6

            if True in ind_sample.values:
                if (np.min(ahe_grain_radius_sample)
                        < ahe_grain_radius_nodes[0, 0]) \
                        or ahe_sample_no == 0:
                    ahe_grain_radius_nodes[:, 0] = \
                        np.min(ahe_grain_radius_sample)
                if (np.max(ahe_grain_radius_sample)
                        > ahe_grain_radius_nodes[0, 1]) \
                        or ahe_sample_no == 0:
                    ahe_grain_radius_nodes[:, 1] = \
                        np.max(ahe_grain_radius_sample)

                # calculate helium production and select min and
                # max values of helium prodcution of all samples
                U = ahe_data['U'][ind_sample].values * 1e-6
                Th = ahe_data['Th'][ind_sample].values * 1e-6
                U238 = (137.88 / 138.88) * U
                U235 = (1.0 / 138.88) * U
                Th232 = Th
                Ur0 = (8 * U238 * decay_constant_238U
                       + 7 * U235 * decay_constant_235U
                       + 6 * Th232 * decay_constant_232Th)

                if np.max(Ur0) > Ur0_max:
                    U_nodes[:, 1] = U[np.argmax(Ur0)]
                    Th_nodes[:, 1] = Th[np.argmax(Ur0)]
                    Ur0_max = Ur0.max()
                if np.max(Ur0) < Ur0_min:
                    U_nodes[:, 0] = U[np.argmin(Ur0)]
                    Th_nodes[:, 0] = Th[np.argmin(Ur0)]
                    Ur0_min = Ur0.min()

        # calculate helium ages for all nodes
        simulated_AHe_data =\
            pybasin_lib.simulate_ahe(
                resample_t, nt_prov, n_nodes, time_array_bp,
                z_nodes, T_nodes, active_nodes,
                prov_start_nodes, prov_end_nodes,
                surface_temp, ahe_grain_radius_nodes, U_nodes, Th_nodes,
                ahe_method=ahe_method,
                alpha=alpha, C0=C0, C1=C1, C2=C2, C3=C3)

        (ahe_age_nodes, ahe_age_nodes_min, ahe_age_nodes_max,
         ahe_node_times_burial, ahe_node_zs) = simulated_AHe_data
    else:
        simulated_AHe_data = None

    nt = T_nodes.shape[0]
    n_ahe_samples = len(ahe_samples_well)

    if n_ahe_samples == 0:
        return (None,
                None,
                None,
                None,
                None,
                simulated_AHe_data)

    # get T history for samples only
    T_ahe_samples = np.zeros((nt, n_ahe_samples))
    for h in range(nt):
        T_ahe_samples[h, :] = \
            np.interp(ahe_samples_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      T_nodes[h, active_nodes[-1]])

    # get burial history of samples
    z_ahe_samples = np.zeros((nt, n_ahe_samples))
    for h in range(nt):
        z_ahe_samples[h, :] = \
            np.interp(ahe_samples_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      z_nodes[h, active_nodes[-1]])

    # get provenance history for samples only
    n_prov = prov_start_nodes.shape[1]
    prov_start_ahe_samples = np.zeros((n_ahe_samples, n_prov))
    prov_end_ahe_samples = np.zeros((n_ahe_samples, n_prov))
    for h in range(n_prov):
        prov_start_ahe_samples[:, h] = \
            np.interp(ahe_samples_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      prov_start_nodes[active_nodes[-1], h])
        prov_end_ahe_samples[:, h] = \
            np.interp(ahe_samples_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      prov_end_nodes[active_nodes[-1], h])

    # get active node array for samples only
    active_nodes_ahe_samples = np.zeros((nt, n_ahe_samples),
                                        dtype=bool)
    for h in range(nt):
        active_nodes_ahe_samples[h, :] = \
            np.interp(ahe_samples_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      active_nodes[h, active_nodes[-1]])

    # assemble grain diameters, U and Th content of each sample
    samples = ahe_samples_well['sample'].values
    ahe_grain_radius_samples = []

    U_samples = []
    Th_samples = []

    for ahe_sample in samples:
        ind_sample = ahe_data['sample'] == ahe_sample
        ahe_grain_radius_sample = ahe_data['grain_radius'][ind_sample].values * 1e-6
        ahe_grain_radius_samples.append(ahe_grain_radius_sample)

        U_sample = ahe_data['U'][ind_sample].values * 1e-6
        U_samples.append(U_sample)
        Th_sample = ahe_data['Th'][ind_sample].values * 1e-6
        Th_samples.append(Th_sample)

    print('calculating AHe for %i samples' % n_ahe_samples)

    simulated_ahe_data_samples =\
        pybasin_lib.simulate_ahe(
            resample_t, nt_prov, n_ahe_samples, time_array_bp,
            z_ahe_samples, T_ahe_samples, active_nodes_ahe_samples,
            prov_start_ahe_samples, prov_end_ahe_samples,
            surface_temp, ahe_grain_radius_samples,
            U_samples, Th_samples,
            ahe_method=ahe_method,
            alpha=alpha, C0=C0, C1=C1, C2=C2, C3=C3)

    (modeled_ahe_age_samples, modeled_ahe_age_samples_min,
     modeled_ahe_age_samples_max, ahe_node_times_burial,
     ahe_node_zs) = simulated_ahe_data_samples

    return (modeled_ahe_age_samples,
            modeled_ahe_age_samples_min,
            modeled_ahe_age_samples_max,
            ahe_node_times_burial,
            ahe_node_zs,
            simulated_AHe_data)


def run_model_and_compare_to_data(well_number, well, well_strat,
                                  strat_info_mod, pybasin_params,
                                  surface_temp, surface_salinity_well,
                                  litho_props,
                                  csv_output_dir,
                                  output_dir,
                                  model_scenario_number,
                                  model_results_series,
                                  T_data, vr_data_df,
                                  aft_samples, aft_ages,
                                  ahe_samples, ahe_data, salinity_data,
                                  vr_method='easyRo',
                                  save_csv_files=True):
    """
    run basin & thermal history model and compare modeled  and observed temperature, salinity,
    vitrinite reflectance, apatite fission track ages and/or apatite (U-Th)/He ages

    """

    # run burial history model
    model_result_vars = \
        pybasin_lib.run_burial_hist_model(well_number, well, well_strat,
                                          strat_info_mod, pybasin_params,
                                          surface_temp, surface_salinity_well,
                                          litho_props,
                                          csv_output_dir,
                                          model_scenario_number,
                                          save_csv_files=save_csv_files)

    if pybasin_params.simulate_salinity is False:
        [geohist_df, time_array, time_array_bp,
         surface_temp_array, basal_hf_array,
         z_nodes, T_nodes, active_nodes,
         n_nodes, n_cells,
         node_strat, node_age,
         prov_start_nodes, prov_end_nodes,porosity_nodes, k_nodes] = \
            model_result_vars
    else:
        [geohist_df, time_array, time_array_bp,
         surface_temp_array, basal_hf_array,
         z_nodes, T_nodes, active_nodes,
         n_nodes, n_cells,
         node_strat, node_age,
         prov_start_nodes, prov_end_nodes,
         C_nodes, surface_salinity_array,
         salinity_lwr_bnd, Dw, q_solute_bottom, q_solute_top] = \
            model_result_vars

    # find out if exhumation end has changed
    exhumed_units = [unit[0] == '-' for unit in geohist_df.index]
    unit_names = list(geohist_df.index)
    exhumed_unit_blocks = []
    start = 0
    while True in exhumed_units[start:]:
        start = start + exhumed_units[start:].index(True)
        end = start + exhumed_units[start:].index(False)
        unit_block = [unit_names[start], unit_names[end-1]]
        exhumed_unit_blocks.append(unit_block)
        start = end

    # calculate cooling during exhumation
    for i, exhumed_unit_block in enumerate(exhumed_unit_blocks):
        start = geohist_df.loc[exhumed_unit_block[-1], 'age_bottom']
        end = geohist_df.loc[exhumed_unit_block[0], 'age_top']

        if end < 0:
            end = 0

        # find temperatures at start and end
        if start * 1e6 in time_array_bp and end * 1e6 in time_array_bp:
            start_ind = np.where(time_array_bp / 1e6 == start)[0][0]
            end_ind = np.where(time_array_bp / 1e6 == end)[0][0]

        else:
            print('could not find exact start and end of '
                  'exhumation in time array')
            print('using closest time instead')
            start_ind = np.argmin(np.abs(time_array_bp / 1e6 - start))
            end_ind = np.argmin(np.abs(time_array_bp / 1e6 - end))

        # calculate cooling
        T_cooling = T_nodes[start_ind] - T_nodes[end_ind]

        # filter out eroded formations
        T_cooling_preserved = T_cooling[active_nodes[end_ind]]

        # store results
        model_results_series['start_exhumation_phase_%i' % i] = start
        model_results_series['end_exhumation_phase_%i' % i] = end
        model_results_series['mean_cooling_exhumation_phase_%i' % i] = T_cooling_preserved.mean()
        model_results_series['min_cooling_exhumation_phase_%i' % i] = T_cooling_preserved.min()
        model_results_series['max_cooling_exhumation_phase_%i' % i] = T_cooling_preserved.max()

    # record max temperature and depth
    model_results_series['max_temperature'] = T_nodes.max()
    model_results_series['max_present_temperature'] = \
        T_nodes[-1, active_nodes[-1]].max()
    model_results_series['max_depth'] = z_nodes.max()

    cebs_input = 'cebs.py'
    if pybasin_params.use_strat_map_input is True \
            and os.path.isfile(cebs_input) is True:

        print('reading model input from cebs.py')

        # from . import cebs
        # model_results_series = cebs.present_temp_in_given_depth(
        #    z_nodes, model_scenario_number, T_nodes, model_results_series)
        # model_results_df = cebs.thermal_conductivity(
        #    model_results_series, node_strat, geohist_df, model_scenario_number,
        #    k_nodes)
        # model_results_df = cebs.porosity(model_results_df, node_strat,
        #                                 geohist_df, model_scenario_number,
        #                                 porosity_nodes)

        msg = "error, the cebs.py module had been deprecated. please contact the developer of the code"
        raise ValueError(msg)

    vr_nodes = None

    ################################
    # simulate vitrinite reflectance
    ################################
    if pybasin_params.simulate_VR is True:

        # find if there are VR samples for this well
        ind = ((vr_data_df['well'] == well)
               & (vr_data_df['depth'] < z_nodes[-1].max()))
        vr_data_well = vr_data_df[ind]

        # interpolate vitrinite reflectance data
        if True in ind.values \
                or pybasin_params.calculate_thermochron_for_all_nodes is True:
            print('calculating vitrinite reflectance for n=%i nodes'
                  % n_nodes)

            vr_nodes = pybasin_lib.calculate_vr(T_nodes,
                                                active_nodes,
                                                time_array,
                                                n_nodes, vr_method=vr_method)

            # store surface and bottom VR value
            model_results_series['vr_surface'] = vr_nodes[-1, active_nodes[-1]][0]
            model_results_series['vr_bottom'] = vr_nodes[-1, active_nodes[-1]][-1]
            
            if pybasin_params.use_strat_map_input is True \
                    and os.path.isfile(cebs_input) is True:
                # model_results_series = cebs.vr_top_bot(
                #      model_results_series, node_strat,
                #     vr_nodes, geohist_df, model_scenario_number)
                # model_results_series = cebs.vr_middle(model_results_series, node_strat,
                #     vr_nodes, geohist_df,
                #     model_scenario_number, z_nodes)
                msg = "error the strat_map_input option has been discontinued"
                raise ValueError(msg)

    if pybasin_params.simulate_salinity is True:
        # store depth to 1 g/L aand 0.035 kg/kg salinity

        # calculate density from salinity and temperature

        density = pybasin_lib.equations_of_state_batzle1992(
            np.zeros_like(T_nodes[-1, active_nodes[-1]]),
            T_nodes[-1, active_nodes[-1]],
            C_nodes[-1, active_nodes[-1]])

        C_nodes_gl = C_nodes[-1, active_nodes[-1]] * density
        depth_to_1g_per_l = np.interp([1.0], C_nodes_gl,
                                      z_nodes[-1, active_nodes[-1]])
        depth_to_seawater_salinity = np.interp([pybasin_params.salinity_seawater],
                                               C_nodes[-1, active_nodes[-1]],
                                               z_nodes[-1, active_nodes[-1]])

        model_results_series['depth_to_C=1gL-1'] = depth_to_1g_per_l
        model_results_series['depth_to_C=0.035kg/kg'] = depth_to_seawater_salinity

    if (pybasin_params.simulate_AFT is True
            or pybasin_params.simulate_AHe is True):
        calculate_thermochron_for_all_nodes = \
            pybasin_params.calculate_thermochron_for_all_nodes
        if pybasin_params.make_model_data_fig is True:
            calculate_thermochron_for_all_nodes = True

    ##############################################################
    # simulate apatite fission track ages and length distributions
    ##############################################################
    simulated_AFT_data = None
    location_has_AFT = False

    if pybasin_params.simulate_AFT is True:

        # find if there is any aft data for this well:
        ind = ((aft_samples['well'] == well) & (aft_samples['depth'] <= z_nodes[-1].max() + 1.0))
        aft_data_well = aft_samples[ind]

        if True in ind.values:
            location_has_AFT = True
        else:
            print('no AFT data found for this location')

        (modeled_aft_age_samples,
         modeled_aft_age_samples_min,
         modeled_aft_age_samples_max,
         aft_ln_mean_samples,
         aft_ln_std_samples,
         aft_sample_times_burial,
         aft_sample_zs,
         aft_sample_times,
         aft_sample_temps,
         simulated_AFT_data,
         z_aft_samples,
         T_samples) = assemble_data_and_simulate_aft(
            pybasin_params.resample_timesteps,
            pybasin_params.provenance_time_nt,
            n_nodes, time_array_bp,
            z_nodes, T_nodes, active_nodes,
            prov_start_nodes, prov_end_nodes,
            pybasin_params.annealing_kinetics_values,
            pybasin_params.annealing_kinetic_param,
            surface_temp,
            aft_data_well,
            calculate_thermochron_for_all_nodes=calculate_thermochron_for_all_nodes,
            annealing_eq=pybasin_params.annealing_equation,
            C0=pybasin_params.C0,
            C1=pybasin_params.C1,
            C2=pybasin_params.C2,
            C3=pybasin_params.C3,
            alpha=pybasin_params.alpha,
            location_has_AFT=location_has_AFT)

        # store surface and bottom VR value
        if simulated_AFT_data is not None:

            (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
             aft_ln_mean_nodes, aft_ln_std_nodes,
             aft_node_times_burial, aft_node_zs,
             aft_node_times, aft_node_temps) = simulated_AFT_data

            model_results_series['aft_age_surface_min'] = aft_age_nodes[active_nodes[-1]][0].min()
            model_results_series['aft_age_surface_max'] = aft_age_nodes[active_nodes[-1]][0].max()
            model_results_series['aft_age_bottom_min'] = aft_age_nodes[active_nodes[-1]][-1].min()
            model_results_series['aft_age_bottom_max'] = aft_age_nodes[ active_nodes[-1]][-1].max()

            calculate_resetting_depth = True
            nodata_val = -99999.9

            # age cutoff for sample/node to be considered fully reset
            # todo: add this as a parameter to pybasin_params.py
            full_resetting_age = 0.5

            if aft_age_nodes_min.min() < full_resetting_age:
                ind_min = aft_age_nodes_min < full_resetting_age
                model_results_series['full_resetting_depth_aft_model_min'] = z_nodes[-1][ind_min][0]
                model_results_series['full_resetting_depth_aft_model_min'] = T_nodes[-1][ind_min][0]
            else:
                model_results_series['full_resetting_depth_aft_model_min'] = nodata_val
                model_results_series['full_resetting_depth_aft_model_min'] = nodata_val

            if aft_age_nodes_max.min() < full_resetting_age:
                ind_max = aft_age_nodes_max < full_resetting_age
                model_results_series['full_resetting_depth_aft_model_max'] = z_nodes[-1][ind_max][0]
                model_results_series['full_resetting_T_aft_model_min'] = T_nodes[-1][ind_max][0]
            else:
                model_results_series['full_resetting_depth_aft_model_max'] = nodata_val
                model_results_series['full_resetting_T_aft_model_min'] = nodata_val

            if calculate_resetting_depth is True and pybasin_params.calculate_thermochron_for_all_nodes is True:

                # modeled resetting depth
                ind_reset_min = aft_age_nodes_min <= node_age
                ind_reset_max = aft_age_nodes_max <= node_age
                if True in ind_reset_min:
                    model_results_series['partial_resetting_depth_aft_model_min'] = \
                        z_nodes[-1][ind_reset_min].min()
                else:
                    model_results_series['partial_resetting_depth_aft_model_min'] = nodata_val
                if True in ind_reset_max:
                    model_results_series['partial_resetting_depth_aft_model_max'] = \
                        z_nodes[-1][ind_reset_max].min()
                else:
                    model_results_series['partial_resetting_depth_aft_model_max'] = nodata_val

    #################################
    # simulate apatite (U-Th)/He ages
    #################################
    location_has_AHe = False

    if pybasin_params.simulate_AHe is True:

        resample_t = pybasin_params.resample_timesteps
        nt_prov = pybasin_params.provenance_time_nt

        # find if there is any aft data for this well:
        location_has_AHe = False
        if ahe_samples is not None:
            ind = ahe_samples['location'] == well
            ahe_samples_well = ahe_samples[ind]

            if True in ind.values:
                location_has_AHe = True

        decay_constant_238U = pybasin_params.decay_constant_238U
        decay_constant_235U = pybasin_params.decay_constant_235U
        decay_constant_232Th = pybasin_params.decay_constant_232Th

        (modeled_ahe_age_samples,
         modeled_ahe_age_samples_min,
         modeled_ahe_age_samples_max,
         ahe_node_times_burial,
         ahe_node_zs,
         simulated_AHe_data) = assemble_data_and_simulate_AHe(
            ahe_samples_well,
            ahe_data,
            decay_constant_238U,
            decay_constant_235U,
            decay_constant_232Th,
            n_nodes,
            resample_t,
            nt_prov,
            time_array_bp,
            z_nodes,
            T_nodes,
            active_nodes,
            prov_start_nodes,
            prov_end_nodes,
            surface_temp,
            calculate_thermochron_for_all_nodes=calculate_thermochron_for_all_nodes,
            C0=pybasin_params.C0,
            C1=pybasin_params.C1,
            C2=pybasin_params.C2,
            C3=pybasin_params.C3,
            alpha=pybasin_params.alpha,
            ahe_method=pybasin_params.ahe_method)

        # store surface and bottom VR value
        if simulated_AHe_data is not None:

            (ahe_age_nodes, ahe_age_nodes_min, ahe_age_nodes_max,
             ahe_node_times_burial, ahe_node_zs) = simulated_AHe_data

            ahe_age_nodes_array = np.array(ahe_age_nodes)
            ahe_age_nodes_min_array = np.min(np.array(ahe_age_nodes_min), axis=1)
            ahe_age_nodes_max_array = np.min(np.array(ahe_age_nodes_min), axis=1)

            model_results_series['ahe_age_surface_min'] = ahe_age_nodes_array[active_nodes[-1]][0].min()
            model_results_series['ahe_age_surface_max'] = ahe_age_nodes_array[active_nodes[-1]][0].max()
            model_results_series['ahe_age_bottom_min'] = ahe_age_nodes_array[active_nodes[-1]][-1].min()
            model_results_series['ahe_age_bottom_max'] = ahe_age_nodes_array[active_nodes[-1]][-1].max()

            calculate_resetting_depth = True
            nodata_val = -99999.9

            # age cutoff for sample/node to be considered fully reset
            # todo: add this as a parameter to pybasin_params.py
            full_resetting_age = 0.5

            if ahe_age_nodes_min_array.min() < full_resetting_age:
                ind_min = ahe_age_nodes_min_array < full_resetting_age
                model_results_series['full_resetting_depth_ahe_model_min'] = z_nodes[-1][ind_min][0]
                model_results_series['full_resetting_T_ahe_model_min'] = T_nodes[-1][ind_min][0]
            else:
                model_results_series['full_resetting_depth_ahe_model_min'] = nodata_val
                model_results_series['full_resetting_T_ahe_model_min'] = nodata_val

            if ahe_age_nodes_max_array.min() < full_resetting_age:
                ind_max = ahe_age_nodes_max_array < full_resetting_age
                model_results_series['full_resetting_depth_ahe_model_max'] = z_nodes[-1][ind_max][0]
                model_results_series['full_resetting_T_ahe_model_max'] = T_nodes[-1][ind_max][0]
            else:
                model_results_series['full_resetting_depth_ahe_model_max'] = nodata_val
                model_results_series['full_resetting_T_ahe_model_max'] = nodata_val

            if calculate_resetting_depth is True and pybasin_params.calculate_thermochron_for_all_nodes is True:

                # modeled resetting depth
                ind_reset_min = ahe_age_nodes_min_array <= node_age
                ind_reset_max = ahe_age_nodes_max_array <= node_age
                if True in ind_reset_min:
                    model_results_series['partial_resetting_depth_ahe_model_min'] = \
                        z_nodes[-1][ind_reset_min].min()
                else:
                    model_results_series['partial_resetting_depth_ahe_model_min'] = nodata_val
                if True in ind_reset_max:
                    model_results_series['partial_resetting_depth_ahe_model_max'] = \
                        z_nodes[-1][ind_reset_max].min()
                else:
                    model_results_series['partial_resetting_depth_ahe_model_max'] = nodata_val

    ##################################
    # calculate model goodness of fit:
    ##################################
    # calculate model error temperature data
    ind = (T_data['well'] == well) & (T_data['depth'] < z_nodes[-1].max())
    T_data_well = T_data[ind]

    if True in ind.values:

        T_gof, T_rmse = model_data_comparison_T(T_data_well, z_nodes,
                                                T_nodes, active_nodes)

        T_model_data = (T_data_well['depth'].values,
                        T_data_well['temperature'].values,
                        T_data_well['temperature_unc_1sigma'].values,
                        T_data_well['data_type'],
                        T_gof, T_rmse)

    else:
        T_rmse = np.nan
        T_gof = np.nan
        T_model_data = None

    # calculate model error VR data
    vr_rmse = np.nan
    vr_gof = np.nan
    if pybasin_params.simulate_VR is True:

        # calculate model error vitrinite reflectance data
        ind = ((vr_data_df['well'] == well)
               & (vr_data_df['depth'] < z_nodes[-1].max()))
        vr_data_well = vr_data_df[ind]

        # interpolate vitrinite reflectance data
        if True in ind.values:

            vr_rmse, vr_gof, vr_data_well = model_data_comparison_VR(
                vr_data_well,
                z_nodes, vr_nodes,
                active_nodes,
                vr_unc_sigma=pybasin_params.vr_unc_sigma)

    # calculate model error AFT data
    aft_age_gof = np.nan
    aft_age_error = np.nan
    if pybasin_params.simulate_AFT is True:

        # calculate model error fission track data
        ind = ((aft_samples['well'] == well) & (aft_samples['depth'] <= z_nodes[-1].max() + 1.0))
        aft_data_well = aft_samples[ind]

        if True in ind.values:

            (aft_age_gof, aft_age_error,
             single_grain_aft_ages,
             single_grain_aft_ages_se_min,
             single_grain_aft_ages_se_plus,
             age_bins, age_pdfs, aft_data_well) = \
                model_data_comparison_AFT_age(aft_data_well, aft_ages,
                                              modeled_aft_age_samples_min,
                                              modeled_aft_age_samples_max)

    # simulate apatite (U-Th)/He data
    ahe_age_gof = np.nan
    ahe_age_error = np.nan

    if pybasin_params.simulate_AHe is True:

        # calculate model error fission track data
        ind = ((ahe_samples['location'] == well) & (ahe_samples['depth'] <= z_nodes[-1].max() + 1.0))
        ahe_samples_well = ahe_samples[ind]

        ahe_age_bin = np.linspace(0, prov_start_nodes.max(), 1000)

        if True in ind.values:

            (ahe_age_gof, ahe_age_error, ahe_ages_all_samples,
             ahe_ages_all_samples_SE,
             ahe_age_bin, ahe_age_pdfs_all_samples, ahe_samples_well) = \
                model_data_comparison_AHe(ahe_samples_well, ahe_data,
                                          ahe_age_bin,
                                          modeled_ahe_age_samples_min,
                                          modeled_ahe_age_samples_max)

    # calculate model error salinity data
    salinity_rmse = np.nan
    salinity_gof = np.nan
    if pybasin_params.simulate_salinity is True:

        ind = (salinity_data['well'] == well) & \
              (salinity_data['depth'] < z_nodes[-1].max())
        salinity_data_well = salinity_data[ind]

        if True in ind.values:
            salinity_gof = model_data_comparison_salinity(
                salinity_data_well, z_nodes, C_nodes, active_nodes)

    # assemble output data
    if pybasin_params.simulate_AFT is True and location_has_AFT is True:

        AFT_data = [simulated_AFT_data,
                    aft_data_well['sample'].values,
                    aft_data_well['depth'].values,
                    aft_data_well['aft_age'].values,
                    aft_data_well['aft_age_stderr_min'].values,
                    aft_data_well['aft_age_stderr_plus'].values,
                    aft_data_well['length_mean'].values,
                    aft_data_well['length_std'].values,
                    modeled_aft_age_samples,
                    single_grain_aft_ages,
                    single_grain_aft_ages_se_min,
                    single_grain_aft_ages_se_plus,
                    age_bins,
                    age_pdfs,
                    aft_age_gof,
                    aft_age_error,
                    aft_sample_times,
                    aft_sample_temps,
                    time_array_bp,
                    z_aft_samples, T_samples,
                    aft_data_well]

    else:
        print('no AFT data found for this location')
        AFT_data = None

    if pybasin_params.simulate_VR is True:
        VR_data = [vr_nodes,
                   vr_data_well['depth'].values,
                   vr_data_well['VR'].values,
                   vr_data_well['VR_min'].values,
                   vr_data_well['VR_max'].values,
                   vr_data_well['VR_unc_1sigma'].values,
                   vr_gof,
                   vr_rmse,
                   vr_data_well]
    else:
        VR_data = None

    if pybasin_params.simulate_salinity is True:
        C_data = [C_nodes, surface_salinity_array,
                  salinity_lwr_bnd,
                  salinity_data_well['depth'].values,
                  salinity_data_well['salinity'].values,
                  salinity_data_well['salinity_unc_1sigma'].values,
                  salinity_rmse,
                  q_solute_bottom, q_solute_top]
    else:
        C_data = None

    if (pybasin_params.simulate_AHe is True
            and location_has_AHe is True):
        AHe_model_data = [ahe_samples_well['depth'].values,
                          ahe_ages_all_samples,
                          ahe_ages_all_samples_SE,
                          ahe_age_bin,
                          ahe_age_pdfs_all_samples,
                          modeled_ahe_age_samples,
                          modeled_ahe_age_samples_min,
                          modeled_ahe_age_samples_max,
                          ahe_age_gof, ahe_age_error,
                          simulated_AHe_data,
                          ahe_samples_well]

    else:
        AHe_model_data = None

    model_run_data = [time_array_bp,
                      surface_temp_array, basal_hf_array,
                      z_nodes, active_nodes, T_nodes,
                      node_strat, node_age]

    return (model_run_data,
            T_model_data, T_gof,
            C_data,
            vr_gof, vr_rmse, VR_data,
            aft_age_gof, aft_age_error, AFT_data,
            ahe_age_gof, ahe_age_error,
            AHe_model_data,
            model_results_series)


def update_model_params_and_run_model_new(model_scenario_number,
                                          pybasin_params,
                                          param_names, param_set,
                                          well_number, well,
                                          model_results_series,
                                          well_strat, well_strat_orig,
                                          strat_info_mod,
                                          surface_temp,
                                          surface_salinity_well,
                                          litho_props,
                                          T_data, vr_data_df,
                                          aft_samples, aft_ages,
                                          ahe_samples, ahe_data,
                                          salinity_data,
                                          csv_output_dir,
                                          output_dir,
                                          log_screen_output,
                                          record_data=True,
                                          save_burial_csv_files=True):

    """
    update the model parameters class and run the model one time

    new version, uses beo style parameter space / sensitivity analysis runs

    returns a tuple with the model result variables
    """

    if log_screen_output is True:
        log_output_dir = os.path.join(output_dir, 'log')
        if os.path.exists(log_output_dir) is False:
            os.mkdir(log_output_dir)

        log_fn = 'log_well_%s_model_scen_%i_PID_%s.out' % (well, model_scenario_number, str(os.getpid()))
        log_path = os.path.join(log_output_dir, log_fn)
        print('redirecting screen output to log file %s' % log_path)

        err_fn = 'error_log_well_%s_model_scen_%i_PID_%s.out' % (well, model_scenario_number, str(os.getpid()))
        err_path = os.path.join(log_output_dir, err_fn)
        print('redirecting screen output for model errors to error log file %s' % err_path)

        sys.stdout = open(log_path, "w")
        sys.stderr = open(err_path, "w")

    well_strat = well_strat_orig.copy()

    # update default parameters in pybasin_params class
    for scenario_param_name, scenario_parameter in zip(param_names, param_set):

        if scenario_parameter is not None:
            # find model parameter name to adjust
            model_param_name = scenario_param_name[:-2]

            if hasattr(pybasin_params, model_param_name) is False:
                msg = 'error, the parameter %s is not in the ModelParameters class in the pybasin_params.py file' \
                      % model_param_name
                msg += ', even though it should be updated for model sensitivity or parameter exploration '
                msg += 'according to the ParameterRanges class. Please check if the spelling of the parameter'

                print(msg)

                raise IndexError(msg)

            print('updating parameter %s from %s to %s'
                  % (model_param_name, str(getattr(pybasin_params, model_param_name)), str(scenario_parameter)))

            # update model parameter
            setattr(pybasin_params, model_param_name, scenario_parameter)

    # check exhumation timing:
    n_exh_phases = len(pybasin_params.exhumation_period_starts)
    pybasin_params.exhumation_period_starts = np.array(pybasin_params.exhumation_period_starts)
    pybasin_params.exhumation_period_ends = np.array(pybasin_params.exhumation_period_ends)

    # set up array for end of exhumation, if not specified directly
    if hasattr(pybasin_params, "exhumation_durations"):
        print('using exhumation duration and using this to calculate end of exhumation phase')
        pybasin_params.exhuamtion_durations = np.array(pybasin_params.exhumation_durations)
        pybasin_params.exhumation_period_ends = pybasin_params.exhumation_period_starts \
                                                - pybasin_params.exhumation_durations
        print('calculated end of exhumation period: ', pybasin_params.exhumation_period_ends)

    ind_nok = pybasin_params.exhumation_period_ends < 0.0
    pybasin_params.exhumation_period_ends[ind_nok] = 0.0

    # calculate end of exhumation
    # for exhumation_phase_id in range(n_exh_phases):
        # exhumation_duration_temp = pybasin_params.exhumation_period_starts[exhumation_phase_id] - \
        #                           pybasin_params.exhumation_period_ends[exhumation_phase_id]

        # exhumation_duration_temp = pybasin_params.exhumation_durations[exhumation_phase_id]

        # check if exhumation duration exceeds starting age of exhumation
        # if (exhumation_duration_temp >= (pybasin_params.exhumation_period_starts[exhumation_phase_id]
        #                                  - pybasin_params.max_hf_timestep) / 1e6):
        #     exhumation_duration_temp = (pybasin_params.exhumation_period_starts[exhumation_phase_id]
        #                                 - (pybasin_params.max_hf_timestep * 3) / 1e6)

        # pybasin_params.exhumation_period_ends[exhumation_phase_id] = \
        #     (pybasin_params.exhumation_period_starts[exhumation_phase_id] - exhumation_duration_temp)

        # print('adjusted duration of exhumation = %0.2f My' \
        #       % exhumation_duration_temp)

    # check if vr method is specified in the parameter file
    if hasattr(pybasin_params, 'vr_method') is False:
        pybasin_params.vr_method = 'easyRo'

    # get values of all input parameters in pybasin_params class
    attributes = inspect.getmembers(
        pybasin_params, lambda attribute: not (inspect.isroutine(attribute)))
    attribute_dict = [attribute for attribute in attributes
                      if not (attribute[0].startswith('__') and
                              attribute[0].endswith('__'))]

    # store input parameter values in dataframe:
    for a in attribute_dict:
        if a[0] in model_results_series.index:
            if type(a[1]) is list or type(a[1]) is np.ndarray:
                model_results_series[a[0]] = str(a[1])
                # store entries of lists or arrays in separate columns
                # TODO: fix this, this does not play nice with the model_results_df in the main loop yet
                # if len(a[1]) < 10:
                #    for i in range(len(a[1])):
                #        col_title = a[0] + '_' + str(i)
                #        model_results_series[col_title] = a[1][i]
            else:
                model_results_series[a[0]] = a[1]

    # run model:
    (model_run_data,
     T_model_data, T_gof,
     C_data,
     vr_gof, vr_rmse, VR_model_data,
     aft_age_gof, aft_age_error, AFT_data,
     ahe_age_gof, ahe_age_error,
     AHe_model_data, model_results_series) = \
        run_model_and_compare_to_data(well_number, well, well_strat,
                                      strat_info_mod, pybasin_params,
                                      surface_temp, surface_salinity_well,
                                      litho_props,
                                      csv_output_dir,
                                      output_dir,
                                      model_scenario_number,
                                      model_results_series,
                                      T_data, vr_data_df,
                                      aft_samples, aft_ages,
                                      ahe_samples, ahe_data,
                                      salinity_data,
                                      vr_method=pybasin_params.vr_method,
                                      save_csv_files=save_burial_csv_files)

    if record_data is True:
        # store gof in model results dataframe
        model_results_series['well'] = well
        model_results_series['T_gof'] = T_gof
        if pybasin_params.simulate_VR is True:
            model_results_series['vr_gof'] = vr_gof
            model_results_series['vr_rmse'] = vr_rmse
        if pybasin_params.simulate_AFT is True:
            model_results_series['aft_age_gof'] = aft_age_gof
            model_results_series['aft_age_error'] = aft_age_error
        if pybasin_params.simulate_AHe is True:
            model_results_series['ahe_gof'] = ahe_age_gof
            model_results_series['ahe_error'] = ahe_age_error

        # calculate cumulative salt loss due to diffusion
        if pybasin_params.simulate_salinity is True:
            (C_nodes, surface_salinity_array, salinity_lwr_bnd,
             salinity_well_depth,
             salinity_well,
             salinity_well_sigma,
             salinity_rmse,
             q_solute_bottom, q_solute_top) = C_data

            (time_array_bp,
             surface_temp_array, basal_hf_array,
             z_nodes, active_nodes, T_nodes,
             node_strat, node_age) = model_run_data

            duration = -np.diff(time_array_bp)
            duration = np.append(duration, time_array_bp[-1])
            duration *= 365.25 * 24 * 60 * 60
            total_solute_flux_top = q_solute_top * duration
            total_solute_flux_bottom = q_solute_bottom * duration

            model_results_series['cumalative_solute_flux_top_kg_m-2'] = np.sum(total_solute_flux_top)
            model_results_series['cumalative_solute_flux_bottom_kg_m-2'] = np.sum(total_solute_flux_bottom)

    # restore well strat file to original
    well_strat = well_strat_orig.copy()

    # screen output GOF data
    print('')
    print('temperature GOF = %0.2f' % T_gof)
    if pybasin_params.simulate_VR is True:
        print('vitrinite reflectance GOF = %0.2f' % vr_gof)
    if pybasin_params.simulate_AFT is True:
        print('AFT age GOF = %0.2f' % aft_age_gof)
        print('AFT age error = %0.2f' % aft_age_error)
    if pybasin_params.simulate_AHe is True:
        print('AHe GOF = %0.2f' % ahe_age_gof)
        print('AHe age error = %0.2f' % ahe_age_error)
    print('')

    return (well_number, well, model_scenario_number,
            model_run_data,
            T_model_data, T_gof,
            C_data,
            vr_gof, VR_model_data,
            aft_age_gof, aft_age_error, AFT_data,
            ahe_age_gof, ahe_age_error,
            AHe_model_data,
            model_results_series)


def check_input_data_files(input_dir, pybasin_params):

    """
    check if all necessary input files are available

    """

    print('checking for input files in %s' % input_dir)

    fns = ['stratigraphy_info.csv', 'well_stratigraphy.csv', 'surface_temperature.csv',
           'lithology_properties.csv', 'temperature_data.csv']

    if pybasin_params.simulate_salinity is True:
        fns += ['surface_salinity.csv', 'salinity_data.csv']

    if pybasin_params.simulate_VR is True:
        fns.append('vitrinite_reflectance.csv')

    if pybasin_params.simulate_AHe is True:
        fns += ['ahe_samples.csv', 'ahe_data.csv']

    if pybasin_params.simulate_AFT is True:
        fns += ['aft_samples.csv', 'aft_data.csv']

    for fn in fns:
        if os.path.exists(os.path.join(input_dir, fn)) is False:
            msg = 'error, could not find input file %s in input directory %s' % (fn, input_dir)
            raise IndexError(msg)

    print('found all necessary input files in %s' % input_dir)

    return


def read_model_input_data(input_dir, pybasin_params):

    """
    Read all model input data .csv files using Pandas

    """

    # read all input data
    print('reading input data')
    # stratigraphy description
    strat_info = pd.read_csv(os.path.join(input_dir, 'stratigraphy_info.csv'), skip_blank_lines=True)
    strat_info = strat_info.set_index('strat_unit')

    # well stratigraphy
    well_strats = pd.read_csv(os.path.join(input_dir, 'well_stratigraphy.csv'), skip_blank_lines=True)

    # surface temperature history
    surface_temp = pd.read_csv(os.path.join(input_dir, 'surface_temperature.csv'), skip_blank_lines=True)

    #
    if pybasin_params.simulate_salinity is True:
        # df_sal = pd.read_csv(os.path.join(input_dir, 'surface_salinity.csv'))
        df_sal_inp = pd.read_csv(os.path.join(input_dir, 'surface_salinity.csv'), skip_blank_lines=True)
        df_sal_inp = df_sal_inp.set_index('well')
        df_sal = df_sal_inp.transpose()
        df_sal['age_start'] = pd.to_numeric(df_sal['age_start'])
        df_sal['age_end'] = pd.to_numeric(df_sal['age_end'])
    else:
        df_sal = None

    # lithology properties
    litho_props = pd.read_csv(os.path.join(input_dir,
                                           'lithology_properties.csv'), skip_blank_lines=True)
    litho_props = litho_props.set_index('lithology')

    # present-day temperature data
    T_data_df = pd.read_csv(os.path.join(input_dir, 'temperature_data.csv'), skip_blank_lines=True)

    if pybasin_params.simulate_VR is True:
        # vitrinite reflectance data
        vr_data_df = pd.read_csv(os.path.join(input_dir,
                                              'vitrinite_reflectance.csv'), skip_blank_lines=True)
    else:
        vr_data_df = None

    # fission track data
    if pybasin_params.simulate_AFT is True:
        aft_samples = pd.read_csv(os.path.join(input_dir, 'aft_samples.csv'), skip_blank_lines=True)

        # load AFT age data
        aft_ages = pd.read_csv(os.path.join(input_dir, 'aft_data.csv'), skip_blank_lines=True)
    else:
        aft_samples = None
        aft_ages = None

    if pybasin_params.simulate_AHe is True:
        ahe_sample_fn = os.path.join(input_dir, 'ahe_samples.csv')
        ahe_data_fn = os.path.join(input_dir, 'ahe_data.csv')
        if os.path.exists(ahe_sample_fn):
            # read apatite U-Th/He (AHe) data
            ahe_samples = pd.read_csv(ahe_sample_fn, skip_blank_lines=True)
        else:
            print('warning, could not find input file %s ' % ahe_sample_fn)
            print('continuing without AHe sample data')
            ahe_samples = None

        if os.path.exists(ahe_data_fn):
            # read apatite U-Th/He (AHe) data
            ahe_data = pd.read_csv(ahe_data_fn, skip_blank_lines=True)
        else:
            print('warning, could not find input file %s ' % ahe_data_fn)
            print('continuing without AHe age data')
            ahe_data = None

    else:
        ahe_samples = None
        ahe_data = None

    if pybasin_params.simulate_salinity is True:
        # read surface salinity bnd condditions
        # Cs = pd.read_csv(os.path.join(input_dir, 'surface_salinity.csv'))

        # and read salinity data
        salinity_data = pd.read_csv(os.path.join(input_dir, 'salinity_data.csv'), skip_blank_lines=True)

    else:
        salinity_data = None

    # T_data, vr_data, aft_samples, aft_ages, ahe_samples, ahe_data, salinity_data

    ########
    # calculate porosity-depth and thermal parameters for each strat unit
    # find lithology columns in stratigraphy dataframe
    cols_temp = strat_info.columns[2:].tolist()
    prov_cols = [col for col in cols_temp if 'provenance' in col]
    litho_cols = [col for col in cols_temp if 'provenance' not in col
                  and 'marine' not in col]

    litho_cols.sort()
    litho_props = litho_props.sort_index()

    # check if lithology data is given for each lithology
    try:
        assert litho_props.index[:-1].tolist() == litho_cols
    except AssertionError as msg:
        print('\nerror, something wrong with input data')
        print('not all lithology units found in strat info file are also in the '
              'lithology_properties file')
        print(msg)
        raise AssertionError(msg)

    # check if no provenance columns left empty
    for p in prov_cols:
        if np.all(strat_info[p].isnull()):
            msg = 'error in parsing stratigraphy_info.csv file, ' \
                  'one or more provenance age columns are empty'
            raise ValueError(msg)

    # check well stratigraphy file
    if well_strats['depth_top'].min() < 0 or well_strats['depth_bottom'].min() < 0:

        msg = 'error, found a negative value for depth in the depth_top or ' \
              'depth_bottom columns in the well stratigraphy file. Please make sure all values for ' \
              'depth are zero or positive'

        raise ValueError(msg)

    # create new copy of dataframe to store results
    strat_info_mod = strat_info.copy()

    # add new lithology properties columns to stratigraphy dataframe
    for litho_prop_name in litho_props.columns:
        strat_info_mod[litho_prop_name] = 0

    # go through all litho properties
    for litho_prop_name in litho_props.columns:

        # go through all lithology types present in strat unit
        for col in litho_cols:
            # find fraction of lithology and multiply by lithology property
            s = strat_info[col].astype(float) * litho_props[litho_prop_name][col]

            # add fraction to column in strat dataframe
            strat_info_mod[litho_prop_name] = strat_info_mod[litho_prop_name] + s

    return (well_strats, strat_info_mod, df_sal,
            T_data_df, vr_data_df,
            aft_samples, aft_ages,
            ahe_samples, ahe_data,
            salinity_data, surface_temp, litho_props)


def select_well_strat(well, well_strats):

    """
    Find the stratigraphy for a particular well

    """

    well_strat = well_strats[well_strats['well'] == well]
    well_strat = well_strat.set_index('strat_unit')

    # copy original well strat file
    well_strat_orig = well_strat.copy()

    if len(well_strat) == 0:
        msg = 'could not find well %s in well strat file' % well
        msg += 'please check your input files and make sure each well can be found in the well stratigraphy file'
        raise IndexError(msg)

    return well_strat, well_strat_orig


def select_well_salinity_bnd(well, salinity_bnd_df):

    """
    Read groundwater salinity values at top boundary
    """

    if well in salinity_bnd_df.columns:
        print('using surface salinity bnd for well %s from file' % well)
        # cols = ['age_start', 'age_end', 'surface_salinity_%s' % well]
        cols = ['age_start', 'age_end', well]
        surface_salinity_well = salinity_bnd_df[cols]
        surface_salinity_well['surface_salinity'] = \
            surface_salinity_well[well]
    else:
        surface_salinity_well = None

    return surface_salinity_well


def setup_model_scenarios_new(ParameterRanges):

    """
    Generate list of model parameters for multiple model runs

    """

    pr = ParameterRanges

    # create list with param values for each model run
    scenario_param_names_raw = dir(pr)
    scenario_param_names = [m for m in scenario_param_names_raw
                            if '__' not in m and '_s' in m]

    scenario_parameter_list = [getattr(pr, p)
                               for p in scenario_param_names]

    # construct list with all parameter combinations
    if pr.parameter_combinations is True:
        scenario_parameter_combinations = \
            list(itertools.product(*scenario_parameter_list))
    else:
        nscens = np.sum(np.array([len(sp) for sp in scenario_parameter_list
                                  if sp is not None]))
        nparams = len(scenario_parameter_list)
        scenario_parameter_combinations = []

        if pr.initial_base_run is True:
            scenario_parameter_combinations.append([None] * nparams)

        for j, sl in enumerate(scenario_parameter_list):
            if sl[0] is not None:
                sc = [None] * nparams
                for sli in sl:
                    sci = list(sc)
                    sci[j] = sli
                    scenario_parameter_combinations.append(sci)

    return scenario_param_names, scenario_parameter_combinations


def setup_model_output_df(n_scenarios):

    """
    Set up a Pandas dataframe to store model results
    """

    model_scenario_numbers = np.arange(n_scenarios)

    cols = ['well',
            'exhumation_magnitude', 'exhumation_start', 'exhumation_duration',
            'exhumation_segment_factor', 'exhumation_duration_factor',
            'basal_heat_flow',
            'T_gof', 'vr_gof', 'aft_age_gof', 'aft_age_error',
            'ahe_gof', 'ahe_error',
            'mean_gof',
            'objective_function',
            'resetting_depth_model_min',
            'resetting_depth_model_max',
            'resetting_depth_data_min',
            'non-resetting_depth_data_max']

    # set up dataframe to store model results
    model_results_df = pd.DataFrame(index=model_scenario_numbers,
                                    columns=cols)
    model_results_df2 = pd.DataFrame(columns=cols, index=[1])

    return model_results_df, model_results_df2


def main():

    parser = argparse.ArgumentParser(description='PyBasin: Model burial, exhumation and thermal history '
                                                 'and low-temperature thermochronology')

    parser.add_argument('model_input_subfolder', metavar='directory', default=None, nargs='?',
                        help='directory containing input dataset for PyBasin')

    parser.add_argument('-w', dest='wells',
                        help='specify wells to include, separated by a comma for multiple wells')

    parser.print_help()

    args = parser.parse_args()

    # check if script dir in python path
    scriptdir = os.path.realpath(sys.path[0])
    if scriptdir not in sys.path:
        sys.path.append(scriptdir)

    if args.model_input_subfolder is not None:
        model_input_subfolder = os.path.join(scriptdir, args.model_input_subfolder)
    else:
        # read default input folder
        fin = open(os.path.join(scriptdir, 'default_input_folder.txt'))
        d = fin.readline()
        fin.close()

        scenario_name = d.split()[-1]
        model_input_subfolder = os.path.join(scriptdir, d.rstrip())

    print('running model input data from folder %s' % model_input_subfolder)

    mpath = os.path.join(model_input_subfolder, 'pybasin_params.py')

    param_module = imp.load_source('pybasin_params', mpath)

    Parameters_original = param_module.ModelParameters
    # model_scenarios = param_module.model_scenarios
    ParameterRanges = param_module.ParameterRanges

    Parameters = Parameters_original()

    year = 365.25 * 24.0 * 60 * 60.

    input_dir = model_input_subfolder
    output_dir = Parameters.output_dir
    datafile_output_dir = Parameters.datafile_output_dir
    csv_output_dir = datafile_output_dir

    if os.path.exists(output_dir) is False:
        os.mkdir(output_dir)

    # pck_output_dir = os.path.join(output_dir, 'model_run_data_files')
    if (Parameters.save_model_run_data is True
            and os.path.exists(datafile_output_dir) is False):
        print('creating directory %s to store model result datafiles'
               % datafile_output_dir)
        os.mkdir(datafile_output_dir)

    fig_output_dir = output_dir
    if os.path.exists(fig_output_dir) is False:
        print('creating directory %s to store model-data comparison figures' \
            % fig_output_dir)
        os.mkdir(fig_output_dir)

    today = datetime.datetime.now()
    today_str = '%i-%i-%i' % (today.day, today.month, today.year)

    check_input_data_files(input_dir, Parameters)

    (well_strats, strat_info_mod, salinity_bnd_df,
     T_data_df, vr_data_df,
     aft_samples, aft_ages,
     ahe_samples, ahe_data,
     salinity_data, surface_temp, litho_props) \
        = read_model_input_data(input_dir, Parameters)

    # model_scenario_param_list, params_to_change = setup_model_scenarios(model_scenarios,
    #                                                  pybasin_params.correct_exhumation_duration)

    model_scenario_param_names, model_scenario_param_list = setup_model_scenarios_new(ParameterRanges)

    if len(model_scenario_param_names) is 0:
        model_scenario_param_names = [None]
        model_scenario_param_list = [[None]]
        print('single model run, setting up parameter set with base case values')

    # check sys arguments to run a particular well
    if args.wells is not None:
        wells_str = args.wells
        wells = wells_str.split(",")
    else:
        wells = Parameters.wells

    print('running the following wells: ', wells)

    n_scenarios = len(wells) * len(model_scenario_param_list)

    # get attributes
    params = inspect.getmembers(
        Parameters, lambda attribute: not (inspect.isroutine(attribute)))
    param_names = [attribute[0] for attribute in params
                       if not (attribute[0].startswith('__') and
                               attribute[0].endswith('__'))]

    # model_results_df, model_results_df2 = setup_model_output_df(n_scenarios)
    columns = ['model_run', 'model_error', 'well', 'computational_time'] + param_names
    columns += ['well',
                'T_gof', 'vr_gof', 'aft_age_gof', 'aft_age_error',
                'ahe_gof', 'ahe_error', 'mean_gof',
                'resetting_depth_model_min',
                'resetting_depth_model_max',
                'resetting_depth_data_min',
                'non-resetting_depth_data_max']

    # set up dataframe to store model results
    index = np.arange(n_scenarios)
    model_results_df = pd.DataFrame(index=index, columns=columns)

    model_scenario_number = 0

    if ParameterRanges.parallel_model_runs is True:
        pool = Pool(processes=ParameterRanges.max_number_of_processes)
        print('initialized parallel model runs with max %i simultaneous processes'
              % ParameterRanges.max_number_of_processes)

        processes = []
        done_processing = []

    start_time = time.time()

    #######################
    # go through all wells
    #######################
    for well_number, well in enumerate(wells):

        print('x' * 20)
        print('well %s, %i/%i' % (well, well_number + 1,
                                  len(wells)))

        if np.any(well_strats['well'] == well) == False:
            raise IOError('error, could not find well %s in well strat file in directory %s' % (well, input_dir))

        well_strat, well_strat_orig = select_well_strat(well, well_strats)

        # read well specific surface salinity bnd condition
        if (Parameters.simulate_salinity is True
                and Parameters.well_specific_surface_salinity_bnd is True):

            surface_salinity_well = select_well_salinity_bnd(well, salinity_bnd_df)
            # else:
            #     print 'could not find well %s in surface salintiy bnd condition file' % well
            #     cols = ['age_start', 'age_end', 'surface_salinity']
            #     surface_salinity_well = salinity_bnd_df[cols]
        else:
            surface_salinity_well = None

        # estimate max number of nodes
        well_strat['thickness'] = well_strat['depth_bottom'] - well_strat['depth_top']
        well_strat['n_nodes_est'] = np.ceil(well_strat['thickness'] / Parameters.max_thickness)
        # buffer of 20000 m for number of nodes to be safe
        n_nodes_est = int(np.sum(well_strat['n_nodes_est']) + 2e4 / Parameters.max_thickness)

        dfc = pd.DataFrame(columns=['well'], index=np.arange(n_nodes_est))
        dfc['well'] = well

        if Parameters.simulate_salinity is True:
            n_ts_est = int(strat_info_mod['age_bottom'].max() * 1e6 / Parameters.max_hf_timestep * 3)
            dfqs = pd.DataFrame(columns=['well'], index=np.arange(n_ts_est))
            dfqs['well'] = well

        save_counter = 0

        # determine interval for saving output to csv file
        save_counter_interval = Parameters.csv_save_interval

        # go through model scenarios specified in the model_scenarios.py file:
        for well_scenario_no, model_scenario_params \
                in enumerate(model_scenario_param_list):

            print('-' * 20)
            print('setting up model scenario %i / %i' % (model_scenario_number + 1, len(model_scenario_param_list)))
            print('-' * 20)

            # restore original parameter values
            Parameters = Parameters_original()

            # estimate total runtime and time left
            if model_scenario_number / 250 == model_scenario_number / 250.0 and model_scenario_number > 0:

                now = time.time()
                time_passed = (now - start_time)
                time_per_scenario = time_passed / model_scenario_number
                time_left = \
                    (n_scenarios - model_scenario_number) * time_per_scenario

                tekst = 'model scenario %i / %i\n' % (model_scenario_number,
                                                      n_scenarios)
                tekst += 'time passed = %s\n' \
                         % datetime.timedelta(seconds=time_passed)
                tekst += 'time left = %s\n' \
                         % datetime.timedelta(seconds=time_left)
                tekst += 'time per scenario = %s\n' \
                         % datetime.timedelta(seconds=time_per_scenario)

                print(tekst)

                print('writing estimated runtime to runtime.txt')

                fout = open('runtime_%s.txt' % well, 'w')
                fout.write(tekst)
                fout.close()

            # restore original well strat dataframe
            well_strat = well_strat_orig.copy()

            model_results_series = model_results_df.loc[model_scenario_number]

            if ParameterRanges.parallel_model_runs is False:

                # run a single model scenario

                log_screen_output = False

                (well_number_check, well_check, model_scenario_number_check,
                 model_run_data,
                 T_model_data, T_gof,
                 C_data,
                 vr_gof, VR_model_data,
                 aft_age_gof, aft_age_error, AFT_data,
                 ahe_age_gof, ahe_age_error,
                 AHe_model_data, model_results_series_updated) = update_model_params_and_run_model_new(
                    model_scenario_number,
                    Parameters,
                    model_scenario_param_names, model_scenario_params,
                    well_number, well,
                    model_results_series,
                    well_strat, well_strat_orig,
                    strat_info_mod,
                    surface_temp,
                    surface_salinity_well,
                    litho_props,
                    T_data_df, vr_data_df,
                    aft_samples, aft_ages,
                    ahe_samples, ahe_data,
                    salinity_data,
                    csv_output_dir,
                    output_dir,
                    log_screen_output)

                well_number_store, well_store, model_scenario_number_store = well_number, well, model_scenario_number

                for col in model_results_series_updated.index:
                    if col not in model_results_df.columns:
                        model_results_df[col] = np.nan

                model_results_df.loc[model_scenario_number_store] = model_results_series_updated

                keep_on_processing = True

            else:

                # set up a new parallel model run
                log_screen_output = True

                p = pool.apply_async(update_model_params_and_run_model_new,
                                     (model_scenario_number,
                                      Parameters,
                                      model_scenario_param_names, model_scenario_params,
                                      well_number, well,
                                      model_results_series,
                                      well_strat, well_strat_orig,
                                      strat_info_mod,
                                      surface_temp,
                                      surface_salinity_well,
                                      litho_props,
                                      T_data_df, vr_data_df,
                                      aft_samples, aft_ages,
                                      ahe_samples, ahe_data,
                                      salinity_data,
                                      csv_output_dir,
                                      output_dir,
                                      log_screen_output))

                processes.append(p)

                done_processing.append(False)

                keep_on_processing = True

            while keep_on_processing is True:
                # if ParameterRanges.parallel_model_runs is False or process_parallel_runs_result is True:

                if ParameterRanges.parallel_model_runs is False:
                    model_result_ready = True

                else:
                    model_result_ready = False

                    # check if there is already a result to process:
                    for ip, p in enumerate(processes):
                        # print('checking if any runs are done')

                        # TODO: this is not so elegant, replace with conditional
                        #  loop instead and wrap output in a seperate function
                        if model_result_ready is False and p.ready() is True and done_processing[ip] is False:

                            print('-' * 20)
                            print('process %i is done' % ip)
                            print('-' * 20)

                            p_result = p.get()

                            (well_number_store, well_store, model_scenario_number_store,
                             model_run_data,
                             T_model_data, T_gof,
                             C_data,
                             vr_gof, VR_model_data,
                             aft_age_gof, aft_age_error, AFT_data,
                             ahe_age_gof, ahe_age_error,
                             AHe_model_data, model_results_series_updated) = p_result

                            for col in model_results_series_updated.index:
                                if col not in model_results_df.columns:
                                    model_results_df[col] = np.nan

                            model_results_df.loc[model_scenario_number_store] = model_results_series_updated

                            model_result_ready = True
                            done_processing[ip] = True

                if model_result_ready is False:
                    pass
                    
                else:
                    # process results of a single model run

                    model_run_data_fig = list(model_run_data)

                    model_run_data_fig.append(T_model_data)
                    model_run_data_fig.append(C_data)
                    model_run_data_fig.append(VR_model_data)
                    model_run_data_fig.append(AFT_data)
                    model_run_data_fig.append(AHe_model_data)

                    today = datetime.datetime.now()
                    today_str = '%i-%i-%i' % (today.day, today.month, today.year)

                    # save salinity and T data
                    (time_array_bp,
                     surface_temp_array, basal_hf_array,
                     z_nodes, active_nodes, T_nodes,
                     node_strat, node_age) = model_run_data

                    # l = len(z_nodes[-1, active_nodes[-1]]) - 1
                    # dfc.loc[:l, 'depth_s%i' % model_scenario_number] = \
                    #    z_nodes[-1, active_nodes[-1]]
                    # dfc['depth_s%i' % model_scenario_number] = z_nodes[-1, active_nodes[-1]]

                    # dfc.loc[:l, 'T_s%i' % model_scenario_number] = \
                    #    T_nodes[-1, active_nodes[-1]]

                    if Parameters.save_model_run_data is True:

                        fn = os.path.join(output_dir,
                                      'model_data_%s_%s_ms%i.pck'
                                      % (well_store, today_str,
                                         model_scenario_number_store))
                        print('saving all data for model run as %s' % fn)
                        fout = open(fn, 'wb')
                        pickle.dump(model_run_data_fig, fout)
                        fout.close()

                    ##############################
                    # save model results .csv file
                    ##############################
                    if save_counter == 0 or save_counter >= save_counter_interval:

                        if wells[0] == wells[-1]:
                            well_txt = wells[0]
                        else:
                            well_txt = '%s-%s' % (wells[0], wells[-1])
                        fn = os.path.join(output_dir, 'model_results_%s_%s_ms0-%i.csv'
                                          % (today_str, well_txt,
                                             n_scenarios))
                        print('saving model results .csv file %s' % fn)
                        model_results_df.to_csv(fn, index_label='model_scenario_number')

                    ####################################
                    # detailed model output to csv files
                    ####################################
                    if Parameters.save_model_run_data is True:

                        # save modeled T-t paths, all nodes
                        _, n_nodes_store = T_nodes.shape

                        rs = Parameters.resample_timesteps
                        n_steps = len(time_array_bp[::rs])

                        cols = ['time_bp']

                        df_tt2 = pd.DataFrame(columns=cols,
                                              index=np.arange(n_steps))

                        df_tt2['time_bp'] = time_array_bp[::rs]

                        for i in range(n_nodes_store):

                            # add depth
                            col_name = 'z_node_%i' % i
                            z_col = z_nodes[::rs, i]
                            z_col[active_nodes[::rs, i] == False] = np.nan
                            df_tt2[col_name] = z_col

                            # add temperature
                            col_name = 'T_node_%i' % i
                            T_col = T_nodes[::rs, i]
                            T_col[active_nodes[::rs, i] == False] = np.nan
                            df_tt2[col_name] = T_col

                        fn = os.path.join(csv_output_dir, 'time_depth_temp_%s_%s_ms%i.csv'
                                          % (well_store, today_str,
                                             model_scenario_number_store))

                        print('saving time-temperature paths to %s' % fn)
                        df_tt2.to_csv(fn, index_label='timestep')

                        # salinity data:
                        if Parameters.simulate_salinity is True:
                            (C_nodes, surface_salinity_array, salinity_lwr_bnd,
                             salinity_well_depth,
                             salinity_well,
                             salinity_well_sigma,
                             salinity_rmse,
                             q_solute_bottom, q_solute_top) = C_data

                            l = len(z_nodes[-1, active_nodes[-1]]) - 1

                            dfc.loc[:l, 'salinity_s%i' % model_scenario_number_store] = \
                                C_nodes[-1, active_nodes[-1]]

                            # save solute flux data
                            columns = ['solute_flux_top_kg_m-2_s-1',
                                      'solute_flux_bottom_kg_m-2_s-1']

                            lqs = len(time_array_bp) - 1

                            dfqs.loc[:lqs, 'solute_flux_top_kg_m-2_s-1_s%i'
                                     % model_scenario_number_store] = q_solute_top
                            dfqs.loc[:lqs, 'solute_flux_bottom_kg_m-2_s-1_s%i'
                                     % model_scenario_number_store] = q_solute_bottom
                            dfqs.loc[:lqs, 'time_yr_s%i' % model_scenario_number_store] = \
                                time_array_bp
                            fn = os.path.join(csv_output_dir,
                                              'solute_flux_data_%s_%s_ms%i.csv'
                                              % (well_store, today_str, n_scenarios))
                            print('saving solute flux data to %s' % fn)
                            dfqs.to_csv(fn, index=False)

                        ## VR and temperature
                        if Parameters.simulate_VR is True:

                            [vr_nodes,
                             vr_depth,
                             vr_obs,
                             vr_min,
                             vr_max,
                             vr_obs_sigma,
                             vr_GOF,
                             vr_rmse,
                             vr_data_well] = VR_model_data

                            if vr_nodes is not None:

                                nn = len(z_nodes[-1, active_nodes[-1]])
                                ind_nn = dfc.index[:nn]

                                dfc.loc[ind_nn, 'depth_s%i' % model_scenario_number_store] = \
                                    z_nodes[-1, active_nodes[-1]]
                                dfc.loc[ind_nn, 'T_s%i' % model_scenario_number_store] = T_nodes[-1, active_nodes[-1]]
                                dfc.loc[ind_nn, 'VR_s%i' % model_scenario_number_store] = vr_nodes[-1, active_nodes[-1]]

                                # save depth vs T and VR data
                                fn = os.path.join(csv_output_dir, 'modeled_depth_T_and_VR_%s_%s_ms%i.csv'
                                                  % (well_store, today_str, model_scenario_number_store))
                                print('saving depth, temperature and VR data to %s' % fn)
                                dfc.to_csv(fn, index=False)

                                # save depth vs T and VR data
                                fn = os.path.join(csv_output_dir, 'model_data_comparison_VR_%s_%s_ms%i.csv'
                                                  % (well_store, today_str, model_scenario_number_store))
                                print('saving depth, temperature and VR data to %s' % fn)
                                vr_data_well.to_csv(fn, index=False)

                        # AFT data:
                        if Parameters.simulate_AFT is True and AFT_data is not None:

                            [simulated_AFT_data,
                             aft_sample,
                             aft_age_depth,
                             aft_age,
                             aft_age_stderr_min,
                             aft_age_stderr_plus,
                             aft_length_mean,
                             aft_length_std,
                             aft_age_samples,
                             single_grain_aft_ages,
                             single_grain_aft_ages_se_min,
                             single_grain_aft_ages_se_plus,
                             aft_age_bins,
                             aft_age_pdfs,
                             aft_age_GOF,
                             aft_age_error,
                             aft_sample_times,
                             aft_sample_temps,
                             time_array_bp,
                             z_aft_samples, T_samples, aft_data_well] = AFT_data

                            # find max number of timesteps for AFT history
                            max_steps = np.max([len(tii) for ti in aft_sample_times for tii in ti])

                            n_models = len(aft_sample_times[0])

                            n_samples = len(aft_sample_times)

                            cols = []
                            for sample_i in range(n_samples):
                                for model_j in range(n_models):
                                    cols = cols + ['name_sample_%i_model_%i' % (sample_i, model_j)]
                                    cols = cols + ['depth_sample_%i_model_%i' % (sample_i, model_j)]
                                    cols = cols + ['aft_age_sample_%i_model_%i' % (sample_i, model_j)]
                                    cols = cols + ['time_sample_%i_model_%i' % (sample_i, model_j)]
                                    cols = cols + ['temp_sample_%i_model_%i' % (sample_i, model_j)]

                            df_tt_aft = pd.DataFrame(columns=cols, index=np.arange(max_steps))

                            for sample_i in range(n_samples):
                                df_tt_aft.loc[0, 'name_sample_%i_model_%i' % (sample_i, model_j)] = \
                                    aft_sample[sample_i]
                                df_tt_aft.loc[0, 'depth_sample_%i_model_%i' % (sample_i, model_j)] = \
                                    aft_age_depth[sample_i]
                                df_tt_aft.loc[0, 'aft_age_sample_%i_model_%i' % (sample_i, model_j)] = \
                                    aft_age[sample_i]

                                for model_j in range(n_models):
                                    steps = len(aft_sample_times[sample_i][model_j])

                                    df_tt_aft.loc[:(steps-1), 'time_sample_%i_model_%i' % (sample_i, model_j)] = \
                                        aft_sample_times[sample_i][model_j]
                                    df_tt_aft.loc[:(steps-1), 'temp_sample_%i_model_%i' % (sample_i, model_j)] = \
                                        aft_sample_temps[sample_i][model_j]

                            # save AFT time-temperature paths
                            fn = os.path.join(csv_output_dir,
                                              'aft_sample_time_temp_%s_%s_ms%i.csv'
                                              % (well_store, today_str,
                                                 model_scenario_number_store))
                            print('saving time-temperature paths AFT samples to %s' % fn)
                            df_tt_aft.to_csv(fn, index=False)

                            ###################################################
                            # save csv file model-data comparison for AFT ages:
                            ###################################################
                            _, n_prov, n_kin = aft_age_samples.shape
                            for sample_i in range(n_samples):
                                for prov_model_i in range(n_prov):
                                    for kin_i in range(n_kin):
                                        aft_data_well.loc[sample_i,
                                                          'modeled_age_prov_%i_kinetic_param_%i'
                                                          % (prov_model_i, kin_i)] = \
                                            aft_age_samples[sample_i, prov_model_i, kin_i]

                            fn = os.path.join(csv_output_dir,
                                              'aft_model_vs_data_%s_%s_ms%i.csv'
                                              % (well_store, today_str,
                                                 model_scenario_number_store))
                            print('saving modeled AFT data for samples to %s' % fn)
                            # df_aft.to_csv(fn)
                            aft_data_well.to_csv(fn)

                        # AHe data
                        if Parameters.simulate_AHe is True and AHe_model_data is not None:

                            [ahe_sample_depths,
                             ahe_ages_all_samples,
                             ahe_ages_all_samples_SE,
                             ahe_age_bin,
                             ahe_age_pdfs,
                             modeled_ahe_age_samples,
                             modeled_ahe_age_samples_min,
                             modeled_ahe_age_samples_max,
                             ahe_age_gof, ahe_age_error,
                             simulated_AHe_data,
                             ahe_data_samples] = AHe_model_data

                            n_samples_ahe = len(modeled_ahe_age_samples)

                            for sample_i, age_mod_sample in enumerate(modeled_ahe_age_samples):
                                for grain_i, age_mod_grains in enumerate(age_mod_sample):
                                    for prov_model_i, age_mod_prov in enumerate(age_mod_grains):
                                        ahe_data_samples.loc[sample_i, 'modeled_ahe_age_grain_%i_prov_%i'
                                                          % (grain_i, prov_model_i)] = age_mod_prov
                                fn = os.path.join(csv_output_dir,
                                                  'ahe_model_vs_data_%s_%s_ms%i.csv'
                                                  % (well_store, today_str,
                                                     model_scenario_number_store))
                                print('saving modeled AHe data to %s' % fn)
                                ahe_data_samples.to_csv(fn)

                    #############################
                    # make a model vs data figure
                    #############################
                    if Parameters.make_model_data_fig is True:
                        fig = pybasin_figures.model_vs_data_figure(
                            model_run_data_fig,
                            contour_variable=Parameters.contour_variable,
                            show_strat_column=Parameters.show_strat_column,
                            show_thermochron_data=Parameters.show_thermochron_data)
                    #    vr_data['depth'], vr_data['VR'], vr_data['unc_range_sigma'])

                        # fn = os.path.join(fig_output_dir,
                        #                   'model_data_fig_%s_%s_ms%i.%s'
                        #                   % (well, today_str,
                        #                      model_scenario_number,
                        #                      pybasin_params.fig_adj))
                        # print 'saving model-data comparison figure %s' % fn
                        # fig.savefig(fn, dpi=200)
                        # pl.clf()

                        if type(Parameters.fig_adj) is list:
                            for fa in Parameters.fig_adj:
                                fn = os.path.join(fig_output_dir,
                                              'model_data_fig_%s_%s_ms%i.%s'
                                              % (well_store, today_str,
                                                 model_scenario_number_store,
                                                 fa))
                                print('saving model-data comparison figure %s' % fn)
                                fig.savefig(fn, dpi=200)

                            pl.clf()
                        else:
                            fn = os.path.join(fig_output_dir,
                                              'model_data_fig_%s_%s_ms%i.%s'
                                              % (well_store, today_str,
                                                 model_scenario_number_store,
                                                 Parameters.fig_adj))
                            print('saving model-data comparison figure %s' % fn)
                            fig.savefig(fn, dpi=200)
                            pl.clf()

                    save_counter += 1

                    if save_counter >= save_counter_interval:
                        save_counter = 0

                if ParameterRanges.parallel_model_runs is True:
                    # n_open_processes = len(processes) - np.sum(np.array([p.ready() for p in processes]))
                    n_open_processes = len(processes) - np.sum(done_processing)
                    last_model_scenario = model_scenario_number >= (n_scenarios - 1)
                    if n_open_processes >= ParameterRanges.max_number_of_processes or last_model_scenario is True:
                        keep_on_processing = True
                    else:
                        keep_on_processing = False

                    if n_open_processes == 0:
                        keep_on_processing = False
                else:
                    keep_on_processing = False

            model_scenario_number += 1
        print('done with all model scenarios')

        # creating separate columns for exhumation rates and heat flow
        cols_to_separate = ['exhumed_thicknesses']
        for col in cols_to_separate:
            # new_df = pd.DataFrame(model_results_df[col].values.tolist(),
            #                      index=model_results_df.index)
            # _, n_new_cols = new_df.shape
            try:
                new_cols = [ast.literal_eval(','.join(c.split())) for c in model_results_df[col].values.tolist()]
                n_new_cols = len(new_cols[0])
                new_cols_a = np.array(new_cols)

                new_col_names = [col + '_' + str(i) for i in range(n_new_cols)]
                for new_col_name, new_col_value in zip(new_col_names, new_cols_a.T):
                    model_results_df[new_col_name] = new_col_value
            except:
                print('failed to separate column %s in model output dataframe / csv file' % col)

        # saving model results to a .csv file
        if wells[0] == wells[-1]:
            well_txt = wells[0]
        else:
            well_txt = '%s-%s' % (wells[0], wells[-1])
        fn = os.path.join(output_dir, 'model_results_%s_%s_ms0-%i_final.csv'
                          % (today_str, well_txt,
                             n_scenarios))
        print('saving model results .csv file %s' % fn)
        model_results_df.to_csv(fn, index_label='model_scenario_number')

    print('done')


if __name__ == '__main__':
    main()
