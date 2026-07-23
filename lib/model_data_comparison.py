"""
functions for comparing modelled and measured data: subsurface
temperature, vitrinite reflectance, apatite fission track and
(U-Th)/He ages, porewater salinity and fluid pressure

each ``model_data_comparison_*`` function takes the observed data for
one data type and the corresponding modelled values, and returns
goodness-of-fit (GOF), root-mean-square error (RMSE) and/or
coefficient of determination (R^2) statistics
"""

import itertools
import logging

import numpy as np
import scipy.stats

from . import pybasin_lib

logger = logging.getLogger(__name__)


def calculate_r_squared(observed, simulated):

    """
    calculate the coefficient of determination (R^2) of a set of
    modeled vs. observed values

    R^2 = 1 - SS_res / SS_tot, with SS_res the sum of squared residuals
    (observed - simulated) and SS_tot the sum of squared deviations of
    the observed values from their mean.

    note this is the plain, unweighted R^2, it does not take into
    account uncertainty in the observed values.
    """

    observed = np.asarray(observed, dtype=float)
    simulated = np.asarray(simulated, dtype=float)

    ind = (np.isnan(observed) == False) & (np.isnan(simulated) == False)
    observed = observed[ind]
    simulated = simulated[ind]

    if len(observed) < 2:
        return np.nan

    ss_res = np.sum((observed - simulated)**2)
    ss_tot = np.sum((observed - np.mean(observed))**2)

    if ss_tot == 0:
        return np.nan

    r_squared = 1.0 - ss_res / ss_tot

    return r_squared


def calculate_r_squared_range(observed, sim_min, sim_max):

    """
    calculate the coefficient of determination (R^2) of a set of
    observed values against a modeled range [sim_min, sim_max], instead
    of a single modeled value

    observed values that fall within the modeled range are assigned an
    error of 0. observed values outside the modeled range are assigned
    an error equal to the distance to the nearest end-member of the
    modeled range.

    note this is the plain, unweighted R^2, it does not take into
    account uncertainty in the observed values.
    """

    observed = np.asarray(observed, dtype=float)
    sim_min = np.asarray(sim_min, dtype=float)
    sim_max = np.asarray(sim_max, dtype=float)

    ind = ((np.isnan(observed) == False)
           & (np.isnan(sim_min) == False)
           & (np.isnan(sim_max) == False))
    observed = observed[ind]
    sim_min = sim_min[ind]
    sim_max = sim_max[ind]

    if len(observed) < 2:
        return np.nan

    # distance to the nearest end-member age, 0 if within the modeled range
    error_below = np.clip(sim_min - observed, 0.0, None)
    error_above = np.clip(observed - sim_max, 0.0, None)
    residual = error_below + error_above

    ss_res = np.sum(residual**2)
    ss_tot = np.sum((observed - np.mean(observed))**2)

    if ss_tot == 0:
        return np.nan

    r_squared = 1.0 - ss_res / ss_tot

    return r_squared


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
    T_data_well.loc[ind_bht_nofit, 'P_fit'] = 0
    T_data_well.loc[ind_bht_nofit, 'residual'] = 15.0
    T_data_well.loc[ind_bht_ok, 'P_fit'] = 1.00
    T_data_well.loc[ind_bht_ok, 'residual'] = 0.0

    T_rmse = np.sqrt(np.mean(T_data_well['residual']**2))
    T_gof = np.mean(T_data_well['P_fit'])
    T_r2 = calculate_r_squared(T_data_well['temperature'],
                               T_data_well['simulated_T'])

    logger.debug('Temperature data:')
    logger.debug(T_data_well)

    return T_gof, T_rmse, T_r2


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
    vr_r2 = calculate_r_squared(vr_data_well['VR'], vr_data_well['simulated_vr'])

    return vr_rmse, vr_gof, vr_r2, vr_data_well


def model_data_comparison_AFT_age(aft_data_well, aft_ages,
                                  modeled_aft_age_samples_min,
                                  modeled_aft_age_samples_max,
                                  verbose=False,
                                  two_sided_gof=False,
                                  gof_age_percentile=(5, 95)):

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

            if two_sided_gof:
                # find user-defined percentile index range of measured age PDF
                pc = np.cumsum(age_pdf)
                pct_low = gof_age_percentile[0] / 100.0
                pct_high = gof_age_percentile[1] / 100.0
                m_start = np.where(pc >= pct_low)[0][0]
                m_end = np.where(pc <= pct_high)[0][-1]
                modeled_n_bins = end_ind - start_ind
                if modeled_n_bins > 0:
                    overlap_bins = max(0, min(end_ind, m_end) - max(start_ind, m_start))
                    reverse_gof = overlap_bins / modeled_n_bins
                    pdf_fit_sum = pdf_fit_sum * reverse_gof

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
            logger.info('no model-data comparison for sample %s, missing age data?' % sample_ix)

    # calculate mean GOF from single grain GOFs for each sample
    aft_age_mean_gof = aft_data_well['GOF_aft_ages'].dropna().mean()
    aft_age_mean_error = aft_data_well['age_error'].dropna().mean()

    # coefficient of determination for modeled vs. measured single grain
    # ages. the model predicts a range of end-member ages per sample
    # (not per grain), so this compares each observed single grain age
    # to the modeled range of the sample it belongs to: grain ages
    # within the modeled range get an error of 0, grain ages outside
    # the range get an error equal to the distance to the nearest
    # end-member age
    single_grain_obs_for_r2 = []
    single_grain_sim_min_for_r2 = []
    single_grain_sim_max_for_r2 = []
    for sample_i, grain_ages in enumerate(single_grain_aft_ages):
        if grain_ages is not None:
            single_grain_obs_for_r2.extend(grain_ages)
            single_grain_sim_min_for_r2.extend(
                [modeled_aft_age_samples_min[sample_i]] * len(grain_ages))
            single_grain_sim_max_for_r2.extend(
                [modeled_aft_age_samples_max[sample_i]] * len(grain_ages))
    aft_age_r2 = calculate_r_squared_range(single_grain_obs_for_r2,
                                           single_grain_sim_min_for_r2,
                                           single_grain_sim_max_for_r2)

    if verbose is True:

        logger.debug(aft_data_well[['sample', 'depth', 'aft_age', 'aft_age_stderr_plus', 'simulated_AFT_min', 'simulated_AFT_max', 'GOF_aft_ages', 'age_error']])

    return (aft_age_mean_gof, aft_age_mean_error, aft_age_r2,
            single_grain_aft_ages, single_grain_aft_ages_se_min,
            single_grain_aft_ages_se_plus,
            age_bins, age_pdfs, aft_data_well)


def model_data_comparison_he(he_samples_well, he_data,
                             he_age_bin,
                             modeled_he_age_samples_min,
                             modeled_he_age_samples_max,
                             two_sided_gof=False,
                             gof_age_percentile=(5, 95)):

    """
    Compare modelled and measure apatite (U-Th)/He ages

    """

    logger.info('calculating GOF He data')

    he_age_pdfs_all_samples = []
    he_ages_all_samples = []
    he_ages_all_samples_SE = []

    # single grain observed vs. modeled ages, for the coefficient of
    # determination. grain ages within the modeled [min, max] range of
    # their grain get an error of 0, grain ages outside this range get
    # an error equal to the distance to the nearest end-member age
    single_grain_obs_he_for_r2 = []
    single_grain_sim_min_he_for_r2 = []
    single_grain_sim_max_he_for_r2 = []

    for he_sample_i, he_sample_ix, he_sample in zip(
            itertools.count(),
            he_samples_well.index,
            he_samples_well['sample']):

        he_age_pdfs = []

        if he_sample in he_data['sample'].values:
            ind_sample = he_data['sample'].values == he_sample

            grain_pdfs = []

            he_ages_all_samples.append(
                he_data.loc[ind_sample, 'he_age_uncorrected'].values)
            he_ages_all_samples_SE.append(
                he_data.loc[ind_sample, 'he_age_uncorrected_se'].values)

            age_error = 0

            for grain_i, he_age_obs, he_age_obs_SE \
                    in zip(itertools.count(),
                           he_data.loc[ind_sample, 'he_age_uncorrected'].values,
                           he_data.loc[ind_sample, 'he_age_uncorrected_se'].values):

                he_age_pdf = scipy.stats.norm.pdf(he_age_bin,
                                                   he_age_obs,
                                                   he_age_obs_SE)

                # normalize to make sum of pdf 1
                he_age_pdf = he_age_pdf / he_age_pdf.sum()

                he_age_pdfs.append(he_age_pdf)

                # find out how much of pdf is covered by simulated
                # end-member U-Th/He ages
                # he_age_pdf = he_age_pdf / he_age_pdf.sum()

                he_age_sim_min = \
                    modeled_he_age_samples_min[he_sample_i][grain_i]
                he_age_sim_max = \
                    modeled_he_age_samples_max[he_sample_i][grain_i]

                single_grain_obs_he_for_r2.append(he_age_obs)
                single_grain_sim_min_he_for_r2.append(he_age_sim_min)
                single_grain_sim_max_he_for_r2.append(he_age_sim_max)

                if he_age_sim_min == 0:
                    start_ind = 0
                else:
                    start_ind = np.where(he_age_sim_min >= he_age_bin)[0][-1]

                if he_age_sim_max == 0.0:
                    end_ind = 0
                else:
                    end_ind = np.where(he_age_sim_max <= he_age_bin)[0][0]

                pdf_fit_sum = np.sum(he_age_pdf[start_ind:end_ind])

                if two_sided_gof:
                    # find user-defined percentile index range of measured He age PDF
                    pc = np.cumsum(he_age_pdf)
                    pct_low = gof_age_percentile[0] / 100.0
                    pct_high = gof_age_percentile[1] / 100.0
                    m_start = np.where(pc >= pct_low)[0][0]
                    m_end = np.where(pc <= pct_high)[0][-1]
                    modeled_n_bins = end_ind - start_ind
                    if modeled_n_bins > 0:
                        overlap_bins = max(0, min(end_ind, m_end) - max(start_ind, m_start))
                        reverse_gof = overlap_bins / modeled_n_bins
                        pdf_fit_sum = pdf_fit_sum * reverse_gof

                grain_pdfs.append(pdf_fit_sum)

                # calculate U-Th/He age error:
                # if np.any(np.isnan(age_pdf)) == False:
                # pc = np.cumsum(he_age_pdf)

                # find +-95% confines of age distribution
                # start_ind = np.where(pc >= 0.05)[0][0]
                # end_ind = np.where(pc <= 0.95)[0][-1]

                age_min = he_age_obs - he_age_obs_SE * 1.96
                age_max = he_age_obs + he_age_obs_SE * 1.96

                # check difference of min modeled aft age and min. value of age distribution
                if he_age_sim_min < age_min:
                    age_error_min = 0
                else:
                    age_error_min = he_age_sim_min - age_min

                # check difference of max modeled aft age and max. value of age distribution
                if he_age_sim_max > age_max:
                    age_error_max = 0
                else:
                    age_error_max = age_max - he_age_sim_max

                age_error += age_error_min + age_error_max
                # aft_data_well.loc[sample_ix, 'age_error'] = age_error

            he_samples_well.loc[he_sample_ix, 'mean_GOF_all_grains'] = \
                np.mean(np.array(grain_pdfs))
            he_samples_well.loc[he_sample_ix, 'min_GOF_all_grains'] = \
                np.min(np.array(grain_pdfs))
            he_samples_well.loc[he_sample_ix, 'max_GOF_all_grains'] = \
                np.max(np.array(grain_pdfs))
            he_samples_well.loc[he_sample_ix, 'mean_he_error'] = age_error / len(grain_pdfs)

        he_age_pdfs_all_samples.append(he_age_pdfs)

    if 'mean_GOF_all_grains' in he_samples_well.columns:
        he_age_gof = he_samples_well['mean_GOF_all_grains'].mean()
    else:
        he_age_gof = np.nan

    if 'mean_he_error' in he_samples_well.columns:
        he_age_error = he_samples_well.loc[he_sample_ix, 'mean_he_error'].mean()
    else:
        he_age_error = 99999.9

    he_age_r2 = calculate_r_squared_range(single_grain_obs_he_for_r2,
                                          single_grain_sim_min_he_for_r2,
                                          single_grain_sim_max_he_for_r2)

    return (he_age_gof, he_age_error, he_age_r2, he_ages_all_samples, he_ages_all_samples_SE,
            he_age_bin, he_age_pdfs_all_samples, he_samples_well)


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
    salinity_r2 = calculate_r_squared(salinity_data_well['salinity'],
                                      salinity_data_well['simulated_salinity'])

    return salinity_gof, salinity_rmse, salinity_r2


def model_data_comparison_pressure(pressure_data_well, z_nodes, P_ex_nodes,
                                   active_nodes, pressure_unc_sigma=1.0):

    """
    Compare modeled and measured fluid (pore) pressure

    modeled pressure is the sum of hydrostatic pressure (1025 kg m-3,
    9.81 m s-2, matching the assumption used throughout the compaction
    and fluid flow model) and modeled excess pressure. measured
    pressure is the formation shut-in pressure (FSIP_MPa column) of a
    drill stem test, compared at the midpoint depth of the tested
    interval (depth_top, depth_bottom columns)

    pressure_unc_sigma is the assumed 1 sigma uncertainty (MPa) of the
    pressure measurements, used unless the input data provides a
    pressure_unc_1sigma column instead
    """

    z_active = z_nodes[-1, active_nodes[-1]]
    hydrostatic_pressure = 1025.0 * 9.81 * z_active / 1.0e6
    total_pressure = hydrostatic_pressure + P_ex_nodes[-1, active_nodes[-1]] / 1.0e6

    pressure_data_well['obs_depth'] = \
        (pressure_data_well['depth_top']
         + pressure_data_well['depth_bottom']) / 2.0

    pressure_data_well['simulated_pressure'] = \
        np.interp(pressure_data_well['obs_depth'], z_active, total_pressure)

    pressure_data_well['residual'] = \
        (pressure_data_well['FSIP_MPa']
         - pressure_data_well['simulated_pressure'])

    if 'pressure_unc_1sigma' in pressure_data_well.columns:
        pressure_sigma = pressure_data_well['pressure_unc_1sigma']
    else:
        pressure_sigma = pressure_unc_sigma

    pressure_data_well['residual_norm'] = \
        pressure_data_well['residual'] / pressure_sigma
    pressure_data_well['P_fit'] = \
        (1.0 - scipy.stats.norm.cdf(
            np.abs(pressure_data_well['residual_norm']))) * 2

    pressure_rmse = np.sqrt(np.mean(pressure_data_well['residual']**2))
    pressure_gof = np.mean(pressure_data_well['P_fit'])
    pressure_r2 = calculate_r_squared(pressure_data_well['FSIP_MPa'],
                                      pressure_data_well['simulated_pressure'])

    return pressure_gof, pressure_rmse, pressure_r2, pressure_data_well
