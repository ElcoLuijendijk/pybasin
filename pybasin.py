"""
PyBasin, version 0.1

recoded in 2014-2015 to simplify the code and to use .csv files
+ pandas for model input

Elco Luijendijk, Goettingen University

<elco.luijendijk@geo.uni-goettingen.de>

"""

import sys
import os
import pdb
import datetime
import pickle
import itertools
import time
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as pl
import scipy.optimize as opt

import lib.pybasin_lib as pybasin_lib
import lib.pybasin_figures as pybasin_figures

# helium diffusion algortihm by Meesters and Dunai (2003)
try:
    import helium_diffusion_models as he
except ImportError:
    print 'warning, failed to import AHe module'


def model_data_comparison_T(T_data_well, z_nodes, T_nodes, active_nodes):

    T_data_well['simulated_T'] = np.interp(T_data_well['depth'],
                                           z_nodes[-1, active_nodes[-1]],
                                           T_nodes[-1, active_nodes[-1]])

    # calculate model error temperature data
    #ind = [T_data_well['temperature_unc_1sigma'].isnull()]
    #T_data_well['temperature_unc_1sigma'][ind] = pybasin_params.vr_unc_sigma
    T_data_well['residual'] = (T_data_well['temperature']
                               - T_data_well['simulated_T'])
    T_data_well['residual_norm'] = (T_data_well['residual']
                                    / T_data_well['temperature_unc_1sigma'])
    T_data_well['P_fit'] = \
        (1.0
         - scipy.stats.norm.cdf(np.abs(T_data_well['residual_norm']))) * 2

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

    return T_gof, T_rmse


def model_data_comparison_VR(vr_data_well, z_nodes, vr_nodes, active_nodes):

    vr_data_well['simulated_vr'] = \
        np.interp(vr_data_well['depth'],
                  z_nodes[-1, active_nodes[-1]],
                  vr_nodes[-1, active_nodes[-1]])

    # calculate model error vitrinite data
    #vr_data_well['VR_unc_1sigma'] = vr_data_well['VR_std']
    ind = vr_data_well['VR_unc_1sigma'].isnull()
    vr_data_well['VR_unc_1sigma'][ind] = pybasin_params.vr_unc_sigma
    vr_data_well['residual'] = (vr_data_well['VR']
                                - vr_data_well['simulated_vr'])
    vr_data_well['residual_norm'] = (vr_data_well['residual']
                                     / vr_data_well['VR_unc_1sigma'])
    vr_data_well['P_fit'] = \
        (1.0
         - scipy.stats.norm.cdf(np.abs(vr_data_well['residual_norm'])))*2

    vr_rmse = np.sqrt(np.mean(vr_data_well['residual']**2))
    vr_gof = np.mean(vr_data_well['P_fit'])

    return vr_rmse, vr_gof


def model_data_comparison_AFT_age(aft_data_well, aft_ages,
                                  modeled_aft_age_samples_min,
                                  modeled_aft_age_samples_max):
    age_bins = []
    age_pdfs = []
    single_grain_aft_ages = []
    single_grain_aft_ages_se_min = []
    single_grain_aft_ages_se_plus = []

    for sample in aft_data_well['sample']:

        if sample in aft_ages['sample'].values:
            ind_sample = aft_ages['sample'].values == sample

            # find single grain ages for this sample
            single_grain_ages_sample = aft_ages['AFT_age'][ind_sample].values
            single_grain_ages_se_min_sample = \
                aft_ages['AFT_age_stderr_min'][ind_sample].values
            single_grain_ages_se_plus_sample = \
                aft_ages['AFT_age_stderr_plus'][ind_sample].values
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
                    aft_data_well['AFT_age'][ind_sample].values,
                    aft_data_well['AFT_age_stderr_min'][ind_sample].values,
                    aft_data_well['AFT_age_stderr_plus'][ind_sample].values)

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
            aft_data_well['simulated_AFT_min'] = modeled_aft_age_samples_min[i]
            aft_data_well['simulated_AFT_max'] = modeled_aft_age_samples_max[i]

            # TODO: find more elegant solution for 0.0 simulated AFT age
            # and check if GOF for AFT ages of 0.0 Ma are correct
            if aft_data_well.ix[sample_ix, 'simulated_AFT_min'] == 0:
                start_ind = 0
            else:
                start_ind = np.where(
                    aft_data_well.ix[sample_ix, 'simulated_AFT_min']
                    >= age_bin)[0][-1]

            if aft_data_well.ix[sample_ix, 'simulated_AFT_max'] == 0.0:
                end_ind = 0
            else:
                # np.where(0.0 >= age_bins)[0]
                end_ind = np.where(
                    aft_data_well.ix[sample_ix, 'simulated_AFT_max']
                    <= age_bin)[0][0]

            pdf_fit_sum = np.sum(age_pdf[start_ind:end_ind])
            pdf_nofit_sum = np.sum(age_pdf[:start_ind]) \
                + np.sum(age_pdf[end_ind:])
            aft_data_well.ix[sample_ix, 'GOF_aft_ages'] = pdf_fit_sum

    # calculate model error:
    for i, sample_ix, age_bin, age_pdf in zip(itertools.count(),
                                              aft_data_well.index,
                                              age_bins,
                                              age_pdfs):

        if np.any(np.isnan(age_pdf)) == False:

            pc = np.cumsum(age_pdf)

            # find +-95% confines of age distribution
            start_ind = np.where(pc >= 0.05)[0][0]
            end_ind = np.where(pc <= 0.95)[0][-1]

            age_min = age_bin[start_ind]
            age_max = age_bin[end_ind]

            # check difference of min modeled aft age and min. value of age distribution
            if modeled_aft_age_samples_min[i] < age_min:
                age_error_min = 0
            else:
                age_error_min = modeled_aft_age_samples_min[i] - age_min

            # check difference of max modeled aft age and max. value of age distribution
            if modeled_aft_age_samples_max[i] > age_max:
                age_error_max = 0
            else:
                age_error_max = age_max - modeled_aft_age_samples_max[i]

            age_error = age_error_min + age_error_max

            aft_data_well.ix[sample_ix, 'age_error'] = age_error

    # calculate mean GOF from single grain GOFs for each sample
    aft_age_gof = aft_data_well['GOF_aft_ages'].dropna().mean()
    aft_age_error = aft_data_well['age_error'].dropna().mean()

    return (aft_age_gof, aft_age_error, single_grain_aft_ages, single_grain_aft_ages_se_min,
            single_grain_aft_ages_se_plus,
            age_bins,
            age_pdfs)


def model_data_comparison_AHe(ahe_samples_well, ahe_data,
                              ahe_age_bin,
                              modeled_ahe_age_samples_min,
                              modeled_ahe_age_samples_max):

    """

    :param ahe_samples_well:
    :param ahe_data:
    :param ahe_age_bin:
    :param modeled_ahe_age_samples_min:
    :param modeled_ahe_age_samples_max:
    :return:
    """

    print 'calculating GOF AHe data'

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
                ahe_data['raw_Ahe_age'][ind_sample].values)
            ahe_ages_all_samples_SE.append(
                ahe_data['raw_Ahe_age_SE'][ind_sample].values)

            age_error = 0

            for grain_i, ahe_age_obs, ahe_age_obs_SE \
                    in zip(itertools.count(),
                           ahe_data['raw_Ahe_age'][ind_sample].values,
                           ahe_data['raw_Ahe_age_SE'][ind_sample].values):

                ahe_age_pdf = scipy.stats.norm.pdf(ahe_age_bin,
                                                   ahe_age_obs,
                                                   ahe_age_obs_SE)

                # normalize to make sum of pdf 1
                ahe_age_pdf = ahe_age_pdf / ahe_age_pdf.sum()

                ahe_age_pdfs.append(ahe_age_pdf)

                # find out how much of pdf is covered by simulated
                # end-member AHe ages
                #ahe_age_pdf = ahe_age_pdf / ahe_age_pdf.sum()

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
                #if np.any(np.isnan(age_pdf)) == False:
                #pdb.set_trace()
                #pc = np.cumsum(ahe_age_pdf)

                # find +-95% confines of age distribution
                #start_ind = np.where(pc >= 0.05)[0][0]
                #end_ind = np.where(pc <= 0.95)[0][-1]

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
                #aft_data_well.ix[sample_ix, 'age_error'] = age_error

        ahe_samples_well.ix[ahe_sample_ix, 'mean_GOF_all_grains'] = np.mean(np.array(grain_pdfs))
        ahe_samples_well.ix[ahe_sample_ix, 'min_GOF_all_grains'] = np.min(np.array(grain_pdfs))
        ahe_samples_well.ix[ahe_sample_ix, 'max_GOF_all_grains'] = np.max(np.array(grain_pdfs))
        ahe_samples_well.ix[ahe_sample_ix, 'ahe_error'] = age_error

        ahe_age_pdfs_all_samples.append(ahe_age_pdfs)

    ahe_age_gof = ahe_samples_well['mean_GOF_all_grains'].mean()
    ahe_age_error = ahe_samples_well.ix[ahe_sample_ix, 'ahe_error'].mean()

    return (ahe_age_gof, ahe_age_error, ahe_ages_all_samples, ahe_ages_all_samples_SE,
            ahe_age_bin, ahe_age_pdfs_all_samples)


def model_data_comparison_salinity(salinity_data_well,
                                   z_nodes, C_nodes, active_nodes):

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
                                   calculate_thermochron_for_all_nodes=False):

    """

    :return:
    """

    #
    #resample_t = pybasin_params.resample_AFT_timesteps
    #nt_prov = pybasin_params.provenance_time_nt
    #pybasin_params.annealing_kinetics_values,
    #pybasin_params.annealing_kinetic_param,

    # pybasin_params.make_model_data_fig is True:


    if calculate_thermochron_for_all_nodes is True:
        # simulate AFT all nodes
        simulated_AFT_data =\
            pybasin_lib.simulate_aft(
                resample_t, nt_prov, n_nodes, time_array_bp,
                z_nodes, T_nodes, active_nodes,
                prov_start_nodes, prov_end_nodes,
                annealing_kinetics_values,
                annealing_kinetic_param,
                surface_temp)

        (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
         aft_ln_mean_nodes, aft_ln_std_nodes,
         aft_node_times_burial, aft_node_zs) = \
            simulated_AFT_data

    else:
        simulated_AFT_data = None

    nt = T_nodes.shape[0]
    n_aft_samples = len(aft_data_well)

    # get T history for samples only
    T_samples = np.zeros((nt, n_aft_samples))
    for h in xrange(nt):
        T_samples[h, :] = \
            np.interp(aft_data_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      T_nodes[h, active_nodes[-1]])

    # get burial history of samples
    z_aft_samples = np.zeros((nt, n_aft_samples))
    for h in xrange(nt):
        z_aft_samples[h, :] = \
            np.interp(aft_data_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      z_nodes[h, active_nodes[-1]])

    # get provenance history for samples only
    n_prov = prov_start_nodes.shape[1]
    prov_start_samples = np.zeros((n_aft_samples, n_prov))
    prov_end_samples = np.zeros((n_aft_samples, n_prov))
    for h in xrange(n_prov):
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
    for h in xrange(nt):
        active_nodes_aft_samples[h, :] = \
            np.interp(aft_data_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      active_nodes[h, active_nodes[-1]])

    # select prov. history for samples:
    simulated_AFT_data_samples =\
        pybasin_lib.simulate_aft(
            resample_t, nt_prov, n_aft_samples, time_array_bp,
            z_aft_samples, T_samples, active_nodes_aft_samples,
            prov_start_samples, prov_end_samples,
            annealing_kinetics_values,
            annealing_kinetic_param,
            surface_temp)

    (modeled_aft_age_samples, modeled_aft_age_samples_min,
     modeled_aft_age_samples_max,
     aft_ln_mean_samples, aft_ln_std_samples,
     aft_sample_times_burial, aft_sample_zs) = simulated_AFT_data_samples

    return (modeled_aft_age_samples,
            modeled_aft_age_samples_min,
            modeled_aft_age_samples_max,
            aft_ln_mean_samples,
            aft_ln_std_samples,
            aft_sample_times_burial,
            aft_sample_zs,
            simulated_AFT_data)


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
                                   calculate_thermochron_for_all_nodes=False):

    """

    :param ahe_data_well:
    :param decay_constant_238U:
    :param decay_constant_235U:
    :param decay_constant_232Th:
    :param n_nodes:
    :param calculate_thermochron_for_all_nodes:
    :return:
    """

    if calculate_thermochron_for_all_nodes is True:

        print '-' * 10
        print 'calculating AHe for all nodes'

        ahe_grain_radius_nodes = np.zeros((n_nodes, 2))
        U_nodes = np.zeros((n_nodes, 2))
        Th_nodes = np.zeros((n_nodes, 2))
        Ur0_max = 0
        Ur0_min = 99999

        # find min and max grain diameters and U and Th contents
        # for this location
        samples = ahe_samples_well['sample'].values
        #ahe_grain_radius_samples = []
        for ahe_sample_no, ahe_sample in enumerate(samples):
            ind_sample = ahe_data['sample'] == ahe_sample
            ahe_grain_radius_sample = \
                ahe_data['grain_radius'][ind_sample].values * 1e-6
            if (np.min(ahe_grain_radius_sample)
                    < ahe_grain_radius_nodes[0, 0]) \
                    or ahe_sample_no==0:
                ahe_grain_radius_nodes[:, 0] = \
                    np.min(ahe_grain_radius_sample)
            if (np.max(ahe_grain_radius_sample)
                    > ahe_grain_radius_nodes[0, 1]) \
                    or ahe_sample_no==0:
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
                surface_temp, ahe_grain_radius_nodes, U_nodes, Th_nodes)

        (ahe_age_nodes, ahe_age_nodes_min, ahe_age_nodes_max,
         ahe_node_times_burial, ahe_node_zs) = simulated_AHe_data
    else:
        simulated_AHe_data = None

    nt = T_nodes.shape[0]
    n_ahe_samples = len(ahe_samples_well)

    # get T history for samples only
    T_ahe_samples = np.zeros((nt, n_ahe_samples))
    for h in xrange(nt):
        T_ahe_samples[h, :] = \
            np.interp(ahe_samples_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      T_nodes[h, active_nodes[-1]])

    # get burial history of samples
    z_ahe_samples = np.zeros((nt, n_ahe_samples))
    for h in xrange(nt):
        z_ahe_samples[h, :] = \
            np.interp(ahe_samples_well['depth'],
                      z_nodes[-1, active_nodes[-1]],
                      z_nodes[h, active_nodes[-1]])

    # get provenance history for samples only
    n_prov = prov_start_nodes.shape[1]
    prov_start_ahe_samples = np.zeros((n_ahe_samples, n_prov))
    prov_end_ahe_samples = np.zeros((n_ahe_samples, n_prov))
    for h in xrange(n_prov):
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
    for h in xrange(nt):
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

    simulated_ahe_data_samples =\
        pybasin_lib.simulate_ahe(
            resample_t, nt_prov, n_ahe_samples, time_array_bp,
            z_ahe_samples, T_ahe_samples, active_nodes_ahe_samples,
            prov_start_ahe_samples, prov_end_ahe_samples,
            surface_temp, ahe_grain_radius_samples,
            U_samples, Th_samples)

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
                                  surface_temp, litho_props,
                                  csv_output_dir,
                                  model_scenario_number,
                                  model_results_df):



    # run burial history model
    model_result_vars = \
        pybasin_lib.run_burial_hist_model(well_number, well, well_strat,
                                          strat_info_mod, pybasin_params,
                                          surface_temp, litho_props,
                                          csv_output_dir,
                                          model_scenario_number)

    if pybasin_params.simulate_salinity is False:
        [geohist_df, time_array, time_array_bp,
         surface_temp_array, basal_hf_array,
         z_nodes, T_nodes, active_nodes,
         n_nodes, n_cells,
         node_strat, node_age,
         prov_start_nodes, prov_end_nodes] = model_result_vars
    else:
        [geohist_df, time_array, time_array_bp,
         surface_temp_array, basal_hf_array,
         z_nodes, T_nodes, active_nodes,
         n_nodes, n_cells,
         node_strat, node_age,
         prov_start_nodes, prov_end_nodes,
         C_nodes, surface_salinity_array,
         salinity_lwr_bnd, Dw] = model_result_vars

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
        start = geohist_df.ix[exhumed_unit_block[-1], 'age_bottom']
        end = geohist_df.ix[exhumed_unit_block[0], 'age_top']

        if end < 0:
            pdb.set_trace()

        # find temperatures at start and end
        if start * 1e6 in time_array_bp and end * 1e6 in time_array_bp:
            start_ind = np.where(time_array_bp / 1e6 == start)[0][0]
            end_ind = np.where(time_array_bp / 1e6 == end)[0][0]

        else:
            print 'warning, could not find exact start and end of ' \
                  'exhumation in time array'
            print 'using closest time instead'
            start_ind = np.argmin(np.abs(time_array_bp / 1e6 - start))
            end_ind = np.argmin(np.abs(time_array_bp / 1e6 - end))

        # calculate cooling
        T_cooling = T_nodes[start_ind] - T_nodes[end_ind]

        # filter out eroded formations
        T_cooling_preserved = T_cooling[active_nodes[end_ind]]

        # store results
        model_results_df.ix[model_scenario_number,
                         'start_exhumation_phase_%i' % i] = start
        model_results_df.ix[model_scenario_number,
                         'end_exhumation_phase_%i' % i] = end
        model_results_df.ix[model_scenario_number,
                         'mean_cooling_exhumation_phase_%i' % i] = \
            T_cooling_preserved.mean()
        model_results_df.ix[model_scenario_number,
                         'min_cooling_exhumation_phase_%i' % i] = \
            T_cooling_preserved.min()
        model_results_df.ix[model_scenario_number,
                         'max_cooling_exhumation_phase_%i' % i] = \
            T_cooling_preserved.max()

    # record max temperature and depth
    model_results_df.ix[model_scenario_number,
                     'max_temperature'] = T_nodes.max()
    model_results_df.ix[model_scenario_number,
                     'max_present_temperature'] = \
        T_nodes[-1, active_nodes[-1]].max()
    model_results_df.ix[model_scenario_number,
                     'max_depth'] = z_nodes.max()

    vr_nodes = None

    ################################
    # simulate vitrinite reflectance
    ################################
    if pybasin_params.simulate_VR is True:

        # find if there are VR samples for this well
        ind = ((vr_data['well'] == well)
               & (vr_data['depth'] < z_nodes[-1].max()))
        vr_data_well = vr_data[ind]

        # interpolate vitrinite reflectance data
        if True in ind.values:
            print 'calculating vitrinite reflectance for n=%i nodes' \
                  % n_nodes

            vr_nodes = pybasin_lib.calculate_vr(T_nodes,
                                                active_nodes,
                                                time_array,
                                                n_nodes)

            # store surface VR value
            model_results_df.ix[model_scenario_number, 'vr_surface'] = \
                vr_nodes[-1, active_nodes[-1]][0]

    #
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
        ind = ((aft_samples['well'] == well)
               & (aft_samples['depth'] <= z_nodes[-1].max() + 1.0))
        aft_data_well = aft_samples[ind]

        if True in ind.values:

            location_has_AFT = True

            (modeled_aft_age_samples,
             modeled_aft_age_samples_min,
             modeled_aft_age_samples_max,
             aft_ln_mean_samples,
             aft_ln_std_samples,
             aft_sample_times_burial,
             aft_sample_zs,
             simulated_AFT_data) = assemble_data_and_simulate_aft(
                pybasin_params.resample_AFT_timesteps,
                pybasin_params.provenance_time_nt,
                n_nodes, time_array_bp,
                z_nodes, T_nodes, active_nodes,
                prov_start_nodes, prov_end_nodes,
                pybasin_params.annealing_kinetics_values,
                pybasin_params.annealing_kinetic_param,
                surface_temp,
                aft_data_well,
                calculate_thermochron_for_all_nodes=
                calculate_thermochron_for_all_nodes)

    #################################
    # simulate apatite (U-Th)/He ages
    #################################
    location_has_AHe = False

    if pybasin_params.simulate_AHe is True:

        resample_t = pybasin_params.resample_AFT_timesteps
        nt_prov = pybasin_params.provenance_time_nt

        # find if there is any aft data for this well:
        ind = ahe_samples['location'] == well
        ahe_samples_well = ahe_samples[ind]

        if True in ind.values:

            location_has_AHe = True

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
                calculate_thermochron_for_all_nodes=
                calculate_thermochron_for_all_nodes)

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
        ind = ((vr_data['well'] == well)
               & (vr_data['depth'] < z_nodes[-1].max()))
        vr_data_well = vr_data[ind]

        # interpolate vitrinite reflectance data
        if True in ind.values:

            vr_rmse, vr_gof = model_data_comparison_VR(vr_data_well,
                                                       z_nodes, vr_nodes,
                                                       active_nodes)

    # calculate model error AFT data
    aft_age_gof = np.nan
    if pybasin_params.simulate_AFT is True:

        # calculate model error fission track data
        ind = ((aft_samples['well'] == well)
           & (aft_samples['depth'] <= z_nodes[-1].max() + 1.0))
        aft_data_well = aft_samples[ind]

        if True in ind.values:
            (aft_age_gof, aft_age_error,
             single_grain_aft_ages,
             single_grain_aft_ages_se_min,
             single_grain_aft_ages_se_plus,
             age_bins, age_pdfs) = \
                model_data_comparison_AFT_age(aft_data_well, aft_ages,
                                              modeled_aft_age_samples_min,
                                              modeled_aft_age_samples_max)

    # simulate apatite (U-Th)/He data
    ahe_age_gof = np.nan
    ahe_age_error = np.nan

    if pybasin_params.simulate_AHe is True:

        # calculate model error fission track data
        ind = ((ahe_samples['location'] == well)
           & (ahe_samples['depth'] <= z_nodes[-1].max() + 1.0))
        ahe_samples_well = ahe_samples[ind]

        ahe_age_bin = np.linspace(0, prov_start_nodes.max(), 100)

        if True in ind.values:

            (ahe_age_gof, ahe_age_error, ahe_ages_all_samples,
             ahe_ages_all_samples_SE,
             ahe_age_bin, ahe_age_pdfs_all_samples) = \
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
    if (pybasin_params.simulate_AFT is True
            and location_has_AFT is True):

        AFT_data = [simulated_AFT_data,
                    aft_data_well['depth'].values,
                    aft_data_well['AFT_age'].values,
                    aft_data_well['AFT_age_stderr_min'].values,
                    aft_data_well['AFT_age_stderr_plus'].values,
                    aft_data_well['length_mean'].values,
                    aft_data_well['length_std'].values,
                    aft_data_well['data_type'].values,
                    modeled_aft_age_samples,
                    single_grain_aft_ages,
                    single_grain_aft_ages_se_min,
                    single_grain_aft_ages_se_plus,
                    age_bins,
                    age_pdfs,
                    aft_age_gof,
                    aft_age_error]

    else:
        AFT_data = None

    if pybasin_params.simulate_VR is True:
        VR_data = [vr_nodes,
                   vr_data_well['depth'].values,
                   vr_data_well['VR'].values,
                   vr_data_well['VR_unc_1sigma'].values,
                   vr_gof]
    else:
        VR_data = None

    if pybasin_params.simulate_salinity is True:
        C_data = [C_nodes, surface_salinity_array,
                  salinity_lwr_bnd,
                  salinity_data_well['depth'].values,
                  salinity_data_well['salinity'].values,
                  salinity_data_well['salinity_unc_1sigma'].values,
                  salinity_rmse]
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
                          simulated_AHe_data]

    else:
        AHe_model_data = None

    model_run_data = [
        time_array_bp,
        surface_temp_array, basal_hf_array,
        z_nodes, active_nodes, T_nodes,
        node_strat, node_age]

    return (model_run_data,
            T_model_data, T_gof,
            C_data,
            vr_gof, VR_data,
            aft_age_gof, aft_age_error, AFT_data,
            ahe_age_gof, ahe_age_error,
            AHe_model_data,
            model_results_df)


def add_one(model_scenario_number):
    return (model_scenario_number + 1)


def update_model_params_and_run_model(model_scenario_params, params_to_change,
                                      model_results_df, model_results_df2,
                                      well_number, well, well_strat,
                                      well_strat_orig,
                                      strat_info_mod, pybasin_params,
                                      surface_temp, litho_props,
                                      csv_output_dir,
                                      output_dir,
                                      model_scenario_number,
                                      return_objective_function,
                                      calibration_target,
                                      record_data,
                                      param_bounds_min,
                                      param_bounds_max):

    well_strat = well_strat_orig.copy()

    (exhumation_magnitude, exhumation_start, exhumation_duration,
     basal_heat_flow, alpha, C1, C2, C3, C4) = (None, None, None, None, None,
                                                None, None, None, None)

    for i in range(len(model_scenario_params)):
        if model_scenario_params[i] > param_bounds_max[i]:
            model_scenario_params[i] = param_bounds_max[i]
        if model_scenario_params[i] < param_bounds_min[i]:
            model_scenario_params[i] = param_bounds_min[i]

    if 'exhumation_magnitude' in params_to_change:
        ind = params_to_change.index('exhumation_magnitude')
        exhumation_magnitude = model_scenario_params[ind]
    if 'exhumation_start' in params_to_change:
        exhumation_start = \
            model_scenario_params[params_to_change.index('exhumation_start')]
    if 'exhumation_duration' in params_to_change:
        exhumation_duration = \
            model_scenario_params[params_to_change.index('exhumation_duration')]
    if 'basal_heat_flow' in params_to_change:
        basal_heat_flow = \
            model_scenario_params[params_to_change.index('basal_heat_flow')]
    if 'AFT_alpha' in params_to_change:
        alpha = model_scenario_params[params_to_change.index('alpha')]
    if 'AFT_C1' in params_to_change:
        C1 = model_scenario_params[params_to_change.index('C1')]
    if 'AFT_C2' in params_to_change:
        C2 = model_scenario_params[params_to_change.index('C2')]
    if 'AFT_C3' in params_to_change:
        C3 = model_scenario_params[params_to_change.index('C3')]
    if 'AFT_C4' in params_to_change:
        C4 = model_scenario_params[params_to_change.index('C4')]

    # update parameter file with new value of exhumation
    exhumation_phase_id = \
        pybasin_params.exhumation_phase_ids.index(
            model_scenarios.exhumation_scenarios_period)
    if exhumation_magnitude is not None:
        pybasin_params.exhumed_thicknesses[exhumation_phase_id] = \
            exhumation_magnitude
        print 'exhumation = %0.1f m' % exhumation_magnitude
    if exhumation_start is not None:
        pybasin_params.exhumation_period_starts[exhumation_phase_id] = \
            exhumation_start
        print 'start of exhumation = %0.2f Ma' % exhumation_start
    if exhumation_duration is not None:
        # check if exhumation duration exceeds starting age of exhumation
        if (exhumation_duration
              > pybasin_params.exhumation_period_starts[exhumation_phase_id]):
            exhumation_duration = \
                (pybasin_params.exhumation_period_starts[exhumation_phase_id]
                 - (pybasin_params.max_hf_timestep * 3) / 1e6)
        pybasin_params.exhumation_period_ends[exhumation_phase_id] = \
            (pybasin_params.exhumation_period_starts[exhumation_phase_id]
             - exhumation_duration)
        print 'duration of exhumation = %0.2f My' % exhumation_duration

    # adjust entire heat flow history
    if basal_heat_flow is not None:
        if model_scenarios.basal_heat_flow_scenario_period == 'all':
            pybasin_params.heatflow_history[:] = basal_heat_flow

            print 'constant basal heat flow = %0.2e' % basal_heat_flow
        elif model_scenarios.basal_heat_flow_scenario_period == 'last':
            pybasin_params.heatflow_history[-1] = basal_heat_flow

            print 'basal heat flow last step = %0.2e' % basal_heat_flow

    if alpha is not None:
        pybasin_params.alpha = alpha
    if C1 is not None:
        pybasin_params.C1 = C1
    if C2 is not None:
        pybasin_params.C2 = C2
    if C3 is not None:
        pybasin_params.C3 = C3
    if C4 is not None:
        pybasin_params.C4 = C4

    if record_data is True:
        # record scenario params in dataframe:
        model_results_df.ix[model_scenario_number, 'exhumation_magnitude'] = \
            pybasin_params.exhumed_thicknesses[exhumation_phase_id]
        model_results_df.ix[model_scenario_number, 'exhumation_start'] = \
            pybasin_params.exhumation_period_starts[exhumation_phase_id]
        model_results_df.ix[model_scenario_number, 'exhumation_duration'] = \
            (pybasin_params.exhumation_period_starts[exhumation_phase_id]
             - pybasin_params.exhumation_period_ends[exhumation_phase_id])
        model_results_df.ix[model_scenario_number, 'basal_heat_flow'] = \
            pybasin_params.heatflow_history[-1]
        model_results_df.ix[model_scenario_number, 'AFT_eq_alpha'] = alpha
        model_results_df.ix[model_scenario_number, 'AFT_eq_C1'] = C1
        model_results_df.ix[model_scenario_number, 'AFT_eq_C2'] = C2
        model_results_df.ix[model_scenario_number, 'AFT_eq_C3'] = C3
        model_results_df.ix[model_scenario_number, 'AFT_eq_C4'] = C4

    (model_run_data,
     T_model_data, T_gof,
     C_data,
     vr_gof, VR_data,
     aft_age_gof, aft_age_error, AFT_data,
     ahe_age_gof, ahe_age_error,
     AHe_model_data, model_results_df) = \
        run_model_and_compare_to_data(well_number, well, well_strat,
                                      strat_info_mod, pybasin_params,
                                      surface_temp, litho_props,
                                      csv_output_dir,
                                      model_scenario_number,
                                      model_results_df)

    if record_data is True:
        # store gof in model results dataframe
        model_results_df.ix[model_scenario_number, 'well'] = well
        model_results_df.ix[model_scenario_number, 'T_gof'] = T_gof
        if pybasin_params.simulate_VR is True:
            model_results_df.ix[model_scenario_number, 'vr_gof'] = vr_gof
        if pybasin_params.simulate_AFT is True:
            model_results_df.ix[model_scenario_number, 'aft_age_gof'] = \
                aft_age_gof
            model_results_df.ix[model_scenario_number, 'aft_age_error'] = \
                aft_age_error
        if pybasin_params.simulate_AHe is True:
            model_results_df.ix[model_scenario_number, 'ahe_gof'] = \
                ahe_age_gof
            model_results_df.ix[model_scenario_number, 'ahe_error'] = \
                ahe_age_error

    if return_objective_function is True:
        # save calibration step input & results:
        pass

    #model_results_df2

    # restore well strat file to original
    well_strat = well_strat_orig.copy()

    # screen output GOF data
    print ''
    print 'temperature GOF = %0.2f' % T_gof
    if pybasin_params.simulate_VR is True:
        print 'vitrinite reflectance GOF = %0.2f' % vr_gof
    if pybasin_params.simulate_AFT is True:
        print 'AFT age GOF = %0.2f' % aft_age_gof
        print 'AFT age error = %0.2f' % aft_age_error
    if pybasin_params.simulate_AHe is True:
        print 'AHe GOF = %0.2f' % ahe_age_gof
        print 'AHe age error = %0.2f' % ahe_age_error
    print ''

    if return_objective_function is True:
        objective_function = 0

        if 'AFT_age' in calibration_target:
            objective_function += aft_age_error
        if 'AHe' in calibration_target:
            objective_function += ahe_age_error

        print 'objective function = %0.2f\n' % objective_function

        #model_results_df2.cols = model_results_df.cols
        ind = model_results_df2.index.max() + 1
        for col in model_results_df.columns:
            model_results_df2.loc[ind, col] = \
                model_results_df.ix[model_scenario_number, col]

        model_results_df2.loc[ind, 'objective_function'] = objective_function

        fnw = '%s_calibration_results.csv' % well
        fn = os.path.join(output_dir, fnw)
        model_results_df2.to_csv(fn)

        return objective_function

    else:

        return (model_run_data,
                T_model_data, T_gof,
                C_data,
                vr_gof, VR_data,
                aft_age_gof, aft_age_error, AFT_data,
                ahe_age_gof, ahe_age_error,
                AHe_model_data)


# check if script dir in python path
scriptdir = os.path.realpath(sys.path[0])
if scriptdir not in sys.path:
    sys.path.append(scriptdir)

print ''

# read default input folder
fin = open(os.path.join(scriptdir, 'default_input_folder.txt'))
d = fin.readline()
fin.close()
scenario_name = d.split()[0]
model_input_subfolder = os.path.join(scriptdir, 'input_data',
                                     scenario_name)
print 'running model input data from folder %s' % model_input_subfolder

# import model parameter and model functions scripts
sys.path.append(model_input_subfolder)
import pybasin_params
import model_scenarios

year = 365.25 * 24.0 * 60 * 60.

decay_constant_238U = 4.916e-18
decay_constant_232Th = 1.57e-18
decay_constant_235U = 3.12e-17

input_dir = pybasin_params.input_dir
output_dir = pybasin_params.output_dir

datafile_output_dir = pybasin_params.datafile_output_dir

if os.path.exists(output_dir) is False:
    os.mkdir(output_dir)

#pck_output_dir = os.path.join(output_dir, 'model_run_data_files')
if (pybasin_params.save_model_run_data is True
    and os.path.exists(datafile_output_dir) is False):
    print 'creating directory %s to store model result datafiles' \
        % datafile_output_dir
    os.mkdir(datafile_output_dir)

csv_output_dir = os.path.join(output_dir, 'burial_history_csv_files')
if os.path.exists(csv_output_dir) is False:
    print 'creating directory %s to store burial history .csv datafiles' \
        % csv_output_dir
    os.mkdir(csv_output_dir)

fig_output_dir = os.path.join(output_dir, 'model_data_fig')
if os.path.exists(fig_output_dir) is False:
    print 'creating directory %s to store model-data comparison figures' \
        % fig_output_dir
    os.mkdir(fig_output_dir)

# read all input data
print 'reading input data'
# stratigraphy description
strat_info = pd.read_csv(os.path.join(input_dir, 'stratigraphy_info.csv'))
strat_info = strat_info.set_index('strat_unit')

# well stratigraphy
well_strats = pd.read_csv(os.path.join(input_dir, 'well_stratigraphy.csv'))

# surface temperature history
surface_temp = pd.read_csv(os.path.join(input_dir, 'surface_temperature.csv'))

# lithology properties
litho_props = pd.read_csv(os.path.join(input_dir, 'lithology_properties.csv'))
litho_props = litho_props.set_index('lithology')

# present-day temperature data
T_data = pd.read_csv(os.path.join(input_dir, 'temperature_data.csv'))

if pybasin_params.simulate_VR is True:
    # vitrinite reflectance data
    vr_data = pd.read_csv(os.path.join(input_dir, 'vitrinite_reflectance.csv'))

# fission track data
if pybasin_params.simulate_AFT is True:
    aft_samples_raw = pd.read_csv(os.path.join(input_dir, 'aft_samples.csv'))

    # filter out samples with low grain count
    ind = ((aft_samples_raw['n_grains'] > pybasin_params.min_grain_no)
            | (aft_samples_raw['n_grains'].isnull()))
    aft_samples = aft_samples_raw[ind]

    # load AFT age data
    aft_ages = pd.read_csv(os.path.join(input_dir, 'aft_ages.csv'))

if pybasin_params.simulate_AHe is True:
    # read apatite U-Th/He (AHe) data
    ahe_samples = pd.read_csv(os.path.join(input_dir, 'ahe_samples.csv'))
    # read apatite U-Th/He (AHe) data
    ahe_data = pd.read_csv(os.path.join(input_dir, 'AHe_data.csv'))

if pybasin_params.simulate_salinity is True:
    # read surface salinity bnd condditions
    #Cs = pd.read_csv(os.path.join(input_dir, 'surface_salinity.csv'))

    # and read salinity data
    salinity_data = pd.read_csv(os.path.join(input_dir, 'salinity_data.csv'))

else:
    Cs = None

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
except AssertionError, msg:
    print '\nerror, something wrong with input data'
    print 'not all lithology units found in strat info file are also in the ' \
          'lithology_properties file'
    print msg
    raise

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

# construct list with exhumation and heatflow estimates
model_scenario_param_list_raw = \
    list(itertools.product(model_scenarios.exhumation_magnitudes,
                           model_scenarios.exhumation_starts,
                           model_scenarios.exhumation_durations,
                           model_scenarios.basal_heat_flow_scenarios))

# go through param list to weed out scenarios with exhumation end < 0 Ma
model_scenario_param_list = []
for i, model_params in enumerate(model_scenario_param_list_raw):
    if model_params[2] <= model_params[1]:
        model_scenario_param_list.append(model_params)
    else:
        pass

# check sys arguments to run a particular well
if len(sys.argv) > 1:
    wells = sys.argv[1:]
else:
    wells = model_scenarios.wells

print 'running the following wells: ', wells

n_scenarios = len(wells) * len(model_scenario_param_list)

model_scenario_numbers = np.arange(n_scenarios)

cols = ['well',
        'exhumation_magnitude', 'exhumation_start', 'exhumation_duration',
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

model_scenario_number = 0

start_time = time.time()

#######################
# go through all wells
#######################
for well_number, well in enumerate(wells):

    print 'x' * 20
    print 'well %s, %i/%i' % (well, well_number + 1,
                              len(wells))

    if np.any(well_strats['well'] == well) is False:
        #print 'warning, well %s not in well strat file'
        raise IOError('error, could not find well %s in well strat file')

    well_strat = well_strats[well_strats['well'] == well]
    well_strat = well_strat.set_index('strat_unit')

    # copy original well strat file
    well_strat_orig = well_strat.copy()


    if pybasin_params.calibrate_model_params is True:

        return_objective_function = True
        record_data = True

        if pybasin_params.load_initial_params is True:

            fn = os.path.join(input_dir, pybasin_params.initial_params_file)
            df_init = pd.read_csv(fn)
            ind = df_init['well'] == well
            df_init_well = df_init.loc[ind]

            # find best model run for current calibration targets
            for i, cal_target in enumerate(pybasin_params.calibration_target):

                if 'AHe' in cal_target:
                    obj_col = 'ahe_gof'
                elif 'AFT' in cal_target:
                    obj_col = 'aft_age_gof'
                elif cal_target == 'T':
                    obj_col = 'T_gof'
                elif cal_target == 'vr':
                    obj_col = 'vr_gof'

                if i == 0:
                    df_init_well['obj_function'] = df_init_well[obj_col]
                else:
                    ind = df_init_well[obj_col] < df_init_well['obj_function']
                    if True in ind.values:
                        df_init_well.loc[ind, 'obj_function'] = \
                            df_init_well.loc[ind, obj_col]

            # copy parameters from best model run
            max_ind = np.argmax(df_init_well['obj_function'])

            model_scenario_params = []

            for param_change in pybasin_params.params_to_change:
                model_scenario_params.append(df_init_well.loc[max_ind, param_change])
        else:
            model_scenario_params = pybasin_params.start_param_values

        bounds = [(minval, maxval) for minval, maxval
                  in zip(pybasin_params.param_bounds_min,
                         pybasin_params.param_bounds_max)]

        args = (pybasin_params.params_to_change, model_results_df,
                model_results_df2,
                well_number, well, well_strat.copy(), well_strat_orig,
                strat_info_mod, pybasin_params,
                surface_temp, litho_props,
                csv_output_dir,
                output_dir,
                model_scenario_number,
                return_objective_function,
                pybasin_params.calibration_target,
                record_data,
                pybasin_params.param_bounds_min,
                pybasin_params.param_bounds_max)

        if (pybasin_params.opt_method == 'L-BFGS-B'
                or pybasin_params.opt_method == 'TNC'
                or pybasin_params.opt_method == 'SLSQP'):

            opt_results = opt.minimize(update_model_params_and_run_model,
                                       model_scenario_params,
                                       args=args,
                                       method=pybasin_params.opt_method,
                                       bounds=bounds)

        else:
             opt_results = opt.minimize(update_model_params_and_run_model,
                                        model_scenario_params,
                                        args=args,
                                        method=pybasin_params.opt_method)

    else:
        # go through all model scenarios:
        for well_scenario_no, model_scenario_params \
                in enumerate(model_scenario_param_list):

            print '-' * 10
            print 'model scenario %i / %i' % (model_scenario_number,
                                              len(model_scenario_param_list))

            # estimate total runtime and time left
            if (model_scenario_number / 250 == model_scenario_number / 250.0
                    and model_scenario_number > 0):

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

                print tekst

                print 'writing estimated runtime to runtime.txt'

                fout = open('runtime_%s.txt' % well, 'w')
                fout.write(tekst)
                fout.close()

            # restore original well strat dataframe
            well_strat = well_strat_orig.copy()
            return_objective_function = False
            record_data = True

            (model_run_data,
                 T_model_data, T_gof,
                 C_data,
                 vr_gof, VR_data,
                 aft_age_gof, aft_age_error, AFT_data,
                 ahe_age_gof, ahe_age_error,
                 AHe_model_data) = update_model_params_and_run_model(
                    model_scenario_params,
                    pybasin_params.params_to_change, model_results_df,
                    model_results_df2,
                    well_number, well, well_strat,
                    strat_info_mod, pybasin_params,
                    surface_temp, litho_props,
                    csv_output_dir,
                    output_dir,
                    model_scenario_number,
                    return_objective_function,
                    pybasin_params.calibration_target,
                    record_data,
                    pybasin_params.param_bounds_min,
                    pybasin_params.param_bounds_max)

            # model_scenario_params, params_to_change,

            # calculate mean goodness of fit of temperature, vitrinite
            # and aft age data
            # TODO: add salinity data...

            # save model run data to .pck file
            #model_run_data = [
            #    time_array_bp,
            #    surface_temp_array, basal_hf_array,
            #    z_nodes, active_nodes, T_nodes,
            #    node_strat, node_age]

            model_run_data_fig = model_run_data

            model_run_data_fig.append(T_model_data)
            model_run_data_fig.append(C_data)
            model_run_data_fig.append(VR_data)
            model_run_data_fig.append(AFT_data)
            model_run_data_fig.append(AHe_model_data)

            today = datetime.datetime.now()
            today_str = '%i-%i-%i' % (today.day, today.month, today.year)

            # save salinity and T data
            if pybasin_params.simulate_salinity is True:

                (time_array_bp,
                 surface_temp_array, basal_hf_array,
                 z_nodes, active_nodes, T_nodes,
                 node_strat, node_age) = model_run_data

                (C_nodes, surface_salinity_array, salinity_lwr_bnd,
                 salinity_well_depth,
                 salinity_well,
                 salinity_well_sigma,
                 salinity_rmse) = C_data

                dfc = pd.DataFrame(columns=['depth', 'salinity', 'T'],
                                   index=np.arange(T_nodes[-1].shape[0]))
                dfc['depth'] = z_nodes[-1]
                dfc['salinity'] = C_nodes[-1]
                dfc['T'] = T_nodes[-1]

                fn = os.path.join(fig_output_dir,
                                  'salinity_data_%s_%s_ms%i.csv'
                                  % (well, today_str, model_scenario_number))

                dfc.to_csv(fn, index=False)

            if pybasin_params.save_model_run_data is True:

                fn = os.path.join(datafile_output_dir,
                                  'model_data_%s_%s_ms%i.pck'
                                  % (well, today_str, model_scenario_number))

                print 'saving model run results to %s' % fn

                fout = open(fn, 'w')
                pickle.dump(model_run_data, fout)
                fout.close()

            #############################
            # make a model vs data figure
            #############################
            if pybasin_params.make_model_data_fig is True:
                fig = pybasin_figures.model_vs_data_figure(
                    model_run_data_fig,
                    contour_variable=pybasin_params.contour_variable)
            #    vr_data['depth'], vr_data['VR'], vr_data['unc_range_sigma'])

                fn = os.path.join(fig_output_dir,
                                  'model_data_fig_%s_%s_ms%i.%s'
                                  % (well, today_str,
                                     model_scenario_number,
                                     pybasin_params.fig_adj))
                print 'saving model-data comparison figure %s' % fn
                fig.savefig(fn, dpi=200)
                pl.clf()

            # save model results .csv file
            if wells[0] == wells[-1]:
                well_txt = wells[0]
            else:
                well_txt = '%s-%s' % (wells[0], wells[-1])
            fn = os.path.join(output_dir, 'model_results_%s_%s_ms0-%i.csv'
                              % (today_str, well_txt,
                                 n_scenarios))
            print 'saving model results .csv file %s' % fn
            model_results_df.to_csv(fn, index_label='model_scenario_number')

            model_scenario_number += 1

print 'done'
