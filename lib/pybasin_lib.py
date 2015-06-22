"""

"""

import pdb
import os
import itertools
import numpy as np
import pandas as pd
import scipy.stats

import easyRo
import AFTlib

try:
    import AFTannealingLib
except ImportError:
    print 'using AFT modules from local lib folder'
    import lib.AFTannealingLib as AFTannealingLib


def integrate_porosity(n0, c, z1, z2):
    
    """
    
    Parameters
    ----------
    n0 : float
        porosity at the surface
    c : float
        compressibility (units ...)
    z1 : float
        depth top of unit (m)
    z2 : float
        depth base of unit (m)
    
    Returns
    -------
    n_int : float
        integrated porosity over depth interval
    """
    
    b_w = n0 / c * (np.exp(-c * z1) - np.exp(-c * z2))
    
    n_int = b_w / (z2 - z1)

    return n_int
    

def calculate_matrix_thickness(n0, c, z1, z2):
    
    """
    Parameters
    ----------
    n0 : float
        porosity at the surface
    c : float
        compressibility (units ...)
    z1 : float
        depth top of unit (m)
    z2 : float
        depth base of unit (m)
    
    Returns
    -------
    b_m : float
        thickness of rock matrix (ie. with pore space removed)
    """
    
    # calculate thickness of pore space bw
    b_w = (z2 - z1) * integrate_porosity(n0, c, z1, z2)
    
    # subtract from total thickness
    b_m = (z2 - z1) - b_w

    return b_m
    

def decompact(bm, n0, c, b0_guess):
    
    """
    
    Parameters
    ----------
    
    b_m : float
        thickness of rock matrix (ie. with pore space removed)
    n0 : float
        porosity at the surface
    c : float
        compressibility (units ...)
    b0_guess : float
        initial guess for decompacted thickness
    
    Returns
    -------
    b0 : float
        new iterative value of decomapcted thickness
    
    """
    
    b0 = bm + (n0 / c) * (1.0 - np.exp(-c * b0_guess))
    
    return b0


def compact(bm, n0, c, z_top, b_guess, max_decompaction_error,
            verbose=False):

    """
    """

    z1 = z_top

    bi = b_guess
    thickness_diff_max = 99999
    n_decompaction_iterations = 0

    while thickness_diff_max > max_decompaction_error:

        z2 = z_top + bi
        bw = integrate_porosity(n0, c, z1, z2) * (z2 - z1)
        bi_new = bw + bm
        bi_diff = bi - bi_new

        thickness_diff_max = np.max(np.abs(bi_diff))

        bi = bi_new

        n_decompaction_iterations += 1

        if verbose is True:
            print 'max change in calculating compacted thickness = %0.2f' \
                % thickness_diff_max

    return bi


def subdivide_strat_units(input_df, max_thickness):

    """
    subdivide strat units:

    """

    # subdivide units that exceed the max thickness
    ind_ths = input_df['present-day_thickness'] > max_thickness

    if np.any(ind_ths.values) == True:

        subdiv_units = []
        subdiv_th = []
        subdiv_age_start = []
        subdiv_age_end = []
        subdiv_bottom = []
        subdiv_top = []

        for ind_th in input_df.index[ind_ths]:

            n_subdiv = int(np.ceil(input_df.ix[ind_th, 'present-day_thickness'] / max_thickness))
            new_th = input_df.ix[ind_th, 'present-day_thickness'] / n_subdiv
            new_duration = (input_df.ix[ind_th, 'age_top'] - input_df.ix[ind_th, 'age_bottom']) / n_subdiv

            age_top = input_df.ix[ind_th, 'age_bottom'] + new_duration
            depth_bottom = input_df.ix[ind_th, 'depth_bottom'] - new_th

            # adjust existing unit
            input_df.ix[ind_th, 'present-day_thickness'] = new_th
            input_df.ix[ind_th, 'age_top'] = input_df.ix[ind_th, 'age_bottom'] + new_duration
            input_df.ix[ind_th, 'depth_top'] = input_df.ix[ind_th, 'depth_bottom'] - new_th

            # add new subdivided units
            for i in xrange(n_subdiv-1):
                subdiv_units.append('%s_s_%i' % (ind_th, i+1))
                subdiv_th.append(new_th)
                subdiv_age_start.append(age_top)
                subdiv_age_end.append(age_top + new_duration)
                subdiv_bottom.append(depth_bottom)
                subdiv_top.append(depth_bottom - new_th)

                age_top += new_duration
                depth_bottom -= new_th

        # add subdivided units to dataframe
        subdiv_df = pd.DataFrame(index=subdiv_units,
                                 columns=input_df.columns,
                                 dtype='float64')

        # copy all values from well strat df
        for i, unit in enumerate(subdiv_units):
            orig_unit = unit.split('_s_')[0]
            try:
                subdiv_df.ix[unit, input_df.columns] = input_df.ix[orig_unit, input_df.columns]
            except KeyError, msg:
                print 'error: ', msg
                pdb.set_trace()

        for i, unit in enumerate(subdiv_units):
            subdiv_df.ix[unit, 'present-day_thickness'] = subdiv_th[i]
            subdiv_df.ix[unit, 'age_bottom'] = subdiv_age_start[i]
            subdiv_df.ix[unit, 'age_top'] = subdiv_age_end[i]
            subdiv_df.ix[unit, 'depth_bottom'] = subdiv_bottom[i]
            subdiv_df.ix[unit, 'depth_top'] = subdiv_top[i]

    # merge original and subdivided dataframes
    output_df = pd.concat((input_df, subdiv_df))

    # sort geohistory to time
    output_df = output_df.sort(['age_bottom'])

    output_df = output_df.convert_objects(convert_numeric=True)

    return output_df


def add_exhumation_phases(well_strat, exhumation_phase_ids,
                          exhumation_period_starts,
                          exhumation_period_ends,
                          exhumed_strat_units,
                          exhumed_thicknesses,
                          original_thicknesses,
                          max_thickness,
                          strat_info_mod):
    
    """
    add exhumation phases to well stratigraphy dataframe
    
    """
    
    well_strat['eroded_thickness'] = 0.0
    #well_strat['later_exhumed_thickness'] = 0.0
    
    new_units_strat_list = []
    new_units_start_list = []
    new_units_end_list = []
    new_units_thickness_list = []
    deposition_codes = []

    oldest_unit = well_strat['age_bottom'].max()

    # add exhumation phases
    for (exhumation_phase_id, exhumation_period_start, exhumation_period_end,
         exhumed_strat_unit, original_thickness, exhumed_thickness) \
            in zip(exhumation_phase_ids,
                   exhumation_period_starts,
                   exhumation_period_ends,
                   exhumed_strat_units,
                   original_thicknesses,
                   exhumed_thicknesses):

        if exhumed_thickness > 0 and exhumation_period_start < oldest_unit:

            wsl = well_strat.index.tolist()

            df_ex = pd.DataFrame(index=exhumed_strat_unit,
                                 columns=['preserved',
                                          'preserved_thickness',
                                          'normal_thickness',
                                          'eroded_thickness'])

            try:
                df_ex['normal_thickness'] = original_thickness
            except ValueError, msg:
                print 'error, original thicknesses list is longer or ' \
                      'shorter than list of strat units in param file'
                print 'check param file'
                raise ValueError(msg)

            df_ex['eroded_thickness'] = 0.0

            # find preserved units and thicknesses
            for unit in df_ex.index:

                unit_preserved1 = [unit in w_unit for w_unit in wsl]
                df_ex.ix[unit, 'preserved'] = np.any(unit_preserved1)
                df_ex.ix[unit, 'preserved_thickness'] = \
                    well_strat['present-day_thickness'][unit_preserved1].sum()

            # list (partly) eroded units
            youngest_unit = df_ex[df_ex['preserved']].index[-1]
            eroded_units = exhumed_strat_unit[exhumed_strat_unit.index(youngest_unit):]

            # calculate eroded thicknesses
            df_ex.ix[eroded_units, 'eroded_thickness'] = \
                df_ex.ix[eroded_units, 'normal_thickness'] - \
                df_ex.ix[eroded_units, 'preserved_thickness']

            # correct < 0 eroded thicknesses
            df_ex['eroded_thickness'][df_ex['eroded_thickness'] < 0] = 0.0

            # calculate total potential exhumation and correct
            if df_ex['eroded_thickness'].sum() < exhumed_thickness:
                df_ex.ix[eroded_units[-1], 'eroded_thickness'] += \
                    exhumed_thickness - df_ex['eroded_thickness'].sum()

            if df_ex['eroded_thickness'].sum() > exhumed_thickness:
                ex_diff = df_ex['eroded_thickness'].sum() - exhumed_thickness
                for eroded_unit in eroded_units[::-1]:
                    if ex_diff > 0:
                        if df_ex.ix[eroded_unit, 'eroded_thickness'] > ex_diff:
                            df_ex.ix[eroded_unit, 'eroded_thickness'] -= ex_diff
                            ex_diff = 0
                        else:
                            ex_diff -= df_ex.ix[eroded_unit, 'eroded_thickness']
                            df_ex.ix[eroded_unit, 'eroded_thickness'] = 0.0

            # recalculate eroded units:
            eroded_units = df_ex.index[df_ex['eroded_thickness']>0].tolist()

            # calculate how many units need to be added:
            df_ex['n_additional_units'] = \
                np.ceil((df_ex['eroded_thickness'] / max_thickness).astype(float))

            df_ex['n_additional_units'][df_ex['n_additional_units'] == 0] = np.nan
            df_ex['thickness_additional_units'] = \
                (df_ex['eroded_thickness'] / df_ex['n_additional_units'])

            # find depositional ages of units
            df_ex.ix[eroded_units, 'age_bottom'] = \
                strat_info_mod.ix[eroded_units, 'age_bottom']
            df_ex.ix[eroded_units, 'age_top'] = \
                strat_info_mod.ix[eroded_units, 'age_top']

            df_ex['duration'] = df_ex['age_bottom'] - df_ex['age_top']

            df_ex.ix[eroded_units, 'duration_preserved_unit'] = \
                ((df_ex.ix[eroded_units, 'preserved_thickness'] /
                  (df_ex.ix[eroded_units, 'preserved_thickness'] +
                   df_ex.ix[eroded_units, 'eroded_thickness']))
                 * df_ex.ix[eroded_units, 'duration'])

            df_ex['duration_non_preserved_units'] = \
                (df_ex['duration'] - df_ex['duration_preserved_unit'])

            df_ex.ix[eroded_units, 'duration_additional_units'] = \
                (df_ex.ix[eroded_units, 'duration_non_preserved_units']
                 / df_ex.ix[eroded_units, 'n_additional_units'])

            # check if end of exhumation not younger than overlying unit
            youngest_unit_index = [i for i, wsli in enumerate(wsl)
                                   if youngest_unit in wsli][0]
            if youngest_unit_index > 0:
                overlying_unit = wsl[youngest_unit_index - 1]
            else:
                overlying_unit = None

            if (overlying_unit is not None
                and exhumation_period_end
                    < well_strat.ix[overlying_unit, 'age_bottom']):
                print 'correcting exhumation end from %0.2f Ma ' \
                      'to age of overlying unit, %0.2f Ma' \
                      % (exhumation_period_end,
                         well_strat.ix[overlying_unit, 'age_bottom'])
                exhumation_period_end = well_strat.ix[overlying_unit,
                                                      'age_bottom']

            if exhumation_period_end < 0:
                print 'correcting exhumation end from %0.2f Ma ' \
                      'to 0 Ma' \
                      % exhumation_period_end
                exhumation_period_end = 0

            # calculate duration of exhumation for each unit
            exhumation_duration_total = (exhumation_period_start
                                         - exhumation_period_end)
            df_ex.ix[eroded_units, 'exhumation_duration_fraction'] = \
                df_ex.ix[eroded_units, 'eroded_thickness'] / exhumed_thickness
            df_ex.ix[eroded_units, 'exhumation_duration'] = \
                (df_ex.ix[eroded_units, 'exhumation_duration_fraction']
                 * exhumation_duration_total)

            df_ex.ix[eroded_units, 'exhumation_units_duration'] = \
                (df_ex.ix[eroded_units, 'exhumation_duration']
                 / df_ex.ix[eroded_units, 'n_additional_units'])

            # distribute exhumation time over all units
            exhumation_start_unit = exhumation_period_start
            for eroded_unit in eroded_units[::-1]:
                df_ex.ix[eroded_unit, 'exhumation_start'] = \
                    exhumation_start_unit
                df_ex.ix[eroded_unit, 'exhumation_end'] = \
                    exhumation_start_unit - \
                    df_ex.ix[eroded_unit, 'exhumation_duration']
                exhumation_start_unit = \
                    exhumation_start_unit - \
                    df_ex.ix[eroded_unit, 'exhumation_duration']

            # adjust deposition time youngest preserved unit
            # find location of youngest preserved unit and subdivisions
            youngest_units = [youngest_unit in w_unit for w_unit in wsl]
            n_youngest_units = np.sum(np.array(youngest_units))

            # update ages youngest preserved units if poartly preserved and
            # partly eroded:
            if df_ex.ix[youngest_unit, 'duration_non_preserved_units'] > 0:
                start_y = well_strat.ix[youngest_units, 'age_bottom'].max()
                #end_y = well_strat.ix[youngest_units, 'age_top'].min()
                durations_fr = -np.arange(n_youngest_units + 1).astype(float) / n_youngest_units
                starts_y = durations_fr[:-1] \
                           * df_ex.ix[youngest_unit, 'duration_preserved_unit'] \
                           + start_y
                ends_y = durations_fr[1:] \
                         * df_ex.ix[youngest_unit, 'duration_preserved_unit'] \
                         + start_y

                well_strat.ix[youngest_units, 'age_bottom'] = starts_y[::-1]
                well_strat.ix[youngest_units, 'age_top'] = ends_y[::-1]

            # create a new dataframe with exhumed units
            for eroded_unit in eroded_units:

                if df_ex.ix[eroded_unit, 'n_additional_units'] > 0:
                    # create lists with timing of deposition new units
                    start = (df_ex.ix[eroded_unit, 'age_bottom'] -
                             df_ex.ix[eroded_unit, 'duration_preserved_unit'])
                    new_units_start = \
                        (-np.arange(
                            df_ex.ix[eroded_unit, 'n_additional_units'])
                         / df_ex.ix[eroded_unit, 'n_additional_units']
                         * df_ex.ix[eroded_unit, 'duration_non_preserved_units']
                         + start)
                    new_units_end = \
                        (new_units_start
                         - df_ex.ix[eroded_unit, 'duration_non_preserved_units']
                         / df_ex.ix[eroded_unit, 'n_additional_units'])

                    # and create list with timing of exhumation of units
                    ex_units_fract = -np.arange(
                            df_ex.ix[eroded_unit, 'n_additional_units'])
                    ex_units_start = \
                        (ex_units_fract
                         * df_ex.ix[eroded_unit, 'exhumation_units_duration']
                         + df_ex.ix[eroded_unit, 'exhumation_start'])[::-1]
                    ex_units_end = \
                        (ex_units_start
                         - df_ex.ix[eroded_unit, 'exhumation_units_duration'])

                    for n_unit, unit_start, unit_end, ex_start, ex_end in zip(
                            itertools.count(),
                            new_units_start,
                            new_units_end,
                            ex_units_start,
                            ex_units_end):

                        # add new units
                        new_units_strat_list.append('+%s_a_%i' % (eroded_unit, n_unit))
                        new_units_start_list.append(unit_start)
                        new_units_end_list.append(unit_end)
                        new_units_thickness_list.append(
                            df_ex.ix[eroded_unit, 'thickness_additional_units'])
                        deposition_codes.append(1)

                        # add exhumation phases
                        new_units_strat_list.append('-%s_a_%i' % (eroded_unit, n_unit))
                        new_units_start_list.append(ex_start)
                        new_units_end_list.append(ex_end)
                        new_units_thickness_list.append(
                            df_ex.ix[eroded_unit, 'thickness_additional_units'])
                        deposition_codes.append(-1)

    # set up new dataframe with new units and exhumation phases
    exhumation_df = pd.DataFrame(index=new_units_strat_list,
                                 columns=well_strat.columns,
                                 dtype='float64')

    exhumation_df['age_bottom'] = np.array(new_units_start_list, dtype=float)
    exhumation_df['age_top'] = np.array(new_units_end_list, dtype=float)
    exhumation_df['deposition_code'] = np.array(deposition_codes,
                                                dtype=float)
    exhumation_df['eroded_thickness'] = np.array(new_units_thickness_list,
                                                 dtype=float)
    exhumation_df['present-day_thickness'] = 0.0

    # merge deposition and exhumation dataframes into one
    output_df = pd.concat((well_strat, exhumation_df))

    # sort new dataframe to have correct sequence:
    output_df = output_df.sort(['age_bottom'])

    return output_df


def find_hiatus(strat_units, age_start, age_end):

    """
    """

    #well_strat_exhumation = pd.DataFrame()
    hiatus_list = []
    hiatus_start_list = []
    hiatus_end_list = []

    for strat_unit_bottom, strat_unit_top, age_end_bottom, age_start_top in \
            zip(strat_units[:-1], strat_units[1:],
                age_end[:-1], age_start[1:]):

        if age_end_bottom > age_start_top + 1e-4:
            print 'found hiatus between units %s (-%0.2f) ' \
                  'and %s (%0.2f-)' \
                  % (strat_unit_bottom,
                     age_end_bottom,
                     strat_unit_top,
                     age_start_top)

            hiatus_list.append('~%s-%s' % (strat_unit_bottom, strat_unit_top))
            hiatus_start_list.append(age_end_bottom)
            hiatus_end_list.append(age_start_top)

        else:
            #print 'no hiatus between units %s (-%0.2f) ' \
            #      'and %s (%0.2f-)' \
            #      % (strat_unit_bottom,
            #         age_end_bottom,
            #         strat_unit_top,
            #         age_start_top)
            pass

    # if last unit deposited before present: add hiatus
    if age_end[-1] > 1e-5:
            hiatus_list.append('~%s-present' % strat_units[-1])
            hiatus_start_list.append(age_end[-1])
            hiatus_end_list.append(0.0)

            print 'added hiatus following depoistion ' \
                  'of youngest unit %s (%0.2f-) ' \
                  % (strat_units[-1], age_end[-1])

    return hiatus_list, hiatus_start_list, hiatus_end_list


def calculate_initial_thickness(max_decompaction_error,
                                surface_porosity, compressibility,
                                depth_top, depth_bottom):

    """

    """

    matrix_thickness = \
        calculate_matrix_thickness(surface_porosity,
                                   compressibility,
                                   depth_top,
                                   depth_bottom)

    initial_thickness = depth_bottom - depth_top

    initial_thickness_diff_max = 99999
    n_decompaction_iterations = 0

    while initial_thickness_diff_max > max_decompaction_error:
        initial_thickness_new = decompact(matrix_thickness, surface_porosity,
                                          compressibility, initial_thickness)
        initial_thickness_diff_max = \
            np.max(np.abs(initial_thickness_new - initial_thickness))
        initial_thickness = initial_thickness_new
        print 'max change in thickness over 1 decompaction iteration %0.2f' \
              % initial_thickness_diff_max
        n_decompaction_iterations += 1

    return initial_thickness


def copy_df_columns(output_df, data_df, rows=None, ignore_age_columns=False):

    """


    Parameters
    ----------

    Returns
    -------

    """

    if rows is None:
        rows = output_df.index

    if ignore_age_columns is False:
        cols = data_df.columns
    else:
        cols_raw = data_df.columns.tolist()
        cols = [col for col in cols_raw
                if ('age_bottom' not in col and 'age_top' not in col)]

    # copy row columns into input dataframe
    for row in rows:
        if row in data_df.index:
            row_name = row
        elif row[0] == '+' and row[1:] in data_df.index:
            print 'adding row info for fully eroded unit %s' % row
            row_name = row[1:]
        elif '_s_' in row and row.split('_s_')[0] in data_df.index:
            print 'adding row info for subdivided unit %s' % row
            row_name = row.split('_s_')[0]
        elif '_a_' in row and row.split('_a_')[0][1:] in data_df.index:
            print 'adding row info for eroded unit %s' % row
            row_name = row.split('_a_')[0][1:]
        else:
            print 'warning, cannot find row item %s in data .csv file' % row
            print 'most likely the names in the well strat file do not ' \
                  'match the names in the stratigraphy data file'
            msg = 'warning, cannot find row item %s in data .csv file. ' \
                  'most likely the names in the well strat file do not ' \
                  'match the names in the stratigraphy data file' % row

            raise KeyError(msg)
            row_name = None

        if row_name is not None:
            output_df.ix[row, cols] = data_df.ix[row_name, cols]

    return output_df


def find_maximum_depth(input_df, exhumation_phases,
                       max_decompaction_error):

    """
    find maximum depth of units by reconstructing burial depth before
    exhumation phases and comparing to present-day depth

    input is a pandas DataFrame with strat units as the index and the
    following columns:
    'depth_top'
    'depth_bottom'
    'eroded_thickness'
    'present-day_thickness'
    'matrix_thickness'
    'surface_porosity'
    'compressibility'

    returns a modified dataframe with the columns 'maximum_depth_top' and
    'maximum_depth_bottom'

    Parameters
    ----------
    input_df : pandas DataFrame object

    exhumation_phases : list
        list with names of exhumation phases
    max_decompaction_error : float


    Returns
    -------

    """

    # add columns for storing max burial depth
    input_df['maximum_depth_top'] = input_df['depth_top']
    input_df['maximum_depth_bottom'] = input_df['depth_bottom']

    input_df['maximum_depth_top_temp'] = np.nan
    input_df['maximum_depth_bottom_temp'] = np.nan
    input_df['maximum_burial_thickness_temp'] = np.nan

    for exhumation_phase in exhumation_phases:

        # find pre-exhumation_strat
        pre_exhumation_start = \
            np.where(input_df.index == exhumation_phase)[0][0]
        pre_exhumation_strat = input_df.index[pre_exhumation_start+1:]

        for i, strat_unit in enumerate(pre_exhumation_strat):

            if strat_unit == pre_exhumation_strat[0]:
                z_top = 0
            else:
                overlying_unit = pre_exhumation_strat[i-1]
                z_top = float(input_df.ix[overlying_unit,
                                          'maximum_depth_bottom_temp'])

            input_df.ix[strat_unit, 'maximum_depth_top_temp'] = z_top

            # add eroded thickness
            if strat_unit[0] == '+':
                input_df.ix[strat_unit, 'maximum_burial_thickness_temp'] = \
                    input_df.ix[strat_unit, 'eroded_thickness']

            # if top is deeper than present-day depth:
            # use present day thickness
            elif input_df.ix[strat_unit, 'maximum_depth_top_temp'] >= \
                    input_df.ix[strat_unit, 'depth_top']:
                input_df.ix[strat_unit, 'maximum_burial_thickness_temp'] = \
                    input_df.ix[strat_unit, 'present-day_thickness']

            # otherwise: calculate compacted thickness:
            else:
                # error in present-day thickness input for SLDNA...
                input_df.ix[strat_unit, 'maximum_burial_thickness_temp'] = \
                    compact(input_df.ix[strat_unit, 'matrix_thickness'],
                            input_df.ix[strat_unit, 'surface_porosity'],
                            input_df.ix[strat_unit, 'compressibility'],
                            z_top,
                            input_df.ix[strat_unit, 'present-day_thickness'],
                            max_decompaction_error)

            input_df.ix[strat_unit, 'maximum_depth_bottom_temp'] = \
                z_top + input_df.ix[strat_unit,
                                    'maximum_burial_thickness_temp']

        # compare to previous estimate of max depth and copy if units more
        # deeply buried
        ind = (input_df['maximum_depth_top'] <
               input_df['maximum_depth_top_temp']) | \
            input_df['maximum_depth_top'].isnull()
        input_df['maximum_depth_top'][ind] = \
            input_df['maximum_depth_top_temp'][ind]
        input_df['maximum_depth_bottom'][ind] = \
            input_df['maximum_depth_bottom_temp'][ind]

    return input_df


def construct_heat_flow_matrix(T, dz, dt, K, rho, c,
                               fixed_upper_temperature,
                               fixed_lower_temperature):

    """
    set up matrix for finite difference eqs.

    """

    nz = len(T)

    # create a diagonal matrix
    A = np.diagflat(np.ones(nz))

    # create the matrix coefficient
    s = K * dt / (rho * c * dz**2)

    # fill matrix
    for i in xrange(1, nz-1):
        A[i, i] = 1 + 2 * s
        A[i, i-1] = -s
        A[i, i+1] = -s

    # upper bnd
    if fixed_upper_temperature is not None:
        A[0, 0] = 1
        A[0, 1] = 0
    else:
        A[0, 0] = 1 + s
        A[0, 1] = -s

    if fixed_lower_temperature is not None:
        A[-1, -1] = 1
        A[-1, -2] = 0
    else:
        #lower bnd
        A[-1, -1] = 1 + s
        A[-1, -2] = -s

    return A


def construct_heat_flow_matrix_variable_z(
        T, z, dt, K, rho, c,
        fixed_upper_temperature,
        fixed_lower_temperature):

    """
    set up matrix for finite difference eqs.

    """

    nz = len(T)

    # create a diagonal matrix
    A = np.diagflat(np.ones(nz))

    # create the matrix coefficient
    s = np.zeros(nz)
    s[1:-1] = 2.0 / (z[2:] - z[:-2])
    s[0] = 1.0 / (z[1] - z[0])
    s[-1] = 1.0 / (z[-1] - z[-2])

    t = np.zeros(nz)
    t[:-1] = K / (z[1:] - z[:-1])

    u = np.zeros(nz)
    u[1:] = K / (z[1:] - z[:-1])

    v = dt / (rho * c)

    # fill matrix
    for i in xrange(1, nz-1):
        A[i, i] = 1 + s[i] * t[i] * v[i] + s[i] * u[i] * v[i]
        A[i, i-1] = -s[i] * u[i] * v[i]
        A[i, i+1] = -s[i] * t[i] * v[i]

    # upper bnd
    if fixed_upper_temperature is not None:
        A[0, 0] = 1
        A[0, 1] = 0
    else:
        A[0, 0] = 1 + s[0] * t[0] * v[0]
        A[0, 1] = -s[0] * t[0] * v[0]

    if fixed_lower_temperature is not None:
        A[-1, -1] = 1
        A[-1, -2] = 0
    else:
        #lower bnd
        A[-1, -1] = 1 + s[-1] * u[-1] * v[-1]
        A[-1, -2] = -s[-1] * u[-1] * v[-1]

    return A


def create_heat_flow_vector(nz, T, Q, dt, rho, c, dz,
                            upper_bnd_flux,
                            lower_bnd_flux,
                            fixed_upper_temperature,
                            fixed_lower_temperature):

    """
    create vector for right hand side implicit FD equation
    """

    # create vector
    b = np.zeros(nz)

    v = dt / (rho * c)

    # fill with values
    b[:] = T[:] + Q * v

    # upper bnd
    if fixed_upper_temperature is not None:
        b[0] = fixed_upper_temperature
    elif upper_bnd_flux is not None:
        #b[0] = T[0] + (Q + upper_bnd_flux / dz) * dt / (rho * c)
        print 'warning upper bnd flux not implemented yet'
        exit()

    # lower bnd
    if fixed_lower_temperature is not None:
        b[-1] = fixed_lower_temperature
    elif lower_bnd_flux is not None:
        s_n = 2.0 / dz
        b[-1] = T[-1] + Q * v + lower_bnd_flux * s_n * v

    return b


def create_heat_flow_vector_variable_z(
        nz, T, Q, dt, rho, c, dz_top, dz_bottom,
        upper_bnd_flux,
        lower_bnd_flux,
        fixed_upper_temperature,
        fixed_lower_temperature):

    """
    create vector for right hand side implicit FD equation
    """

    # create vector
    b = np.zeros(nz)

    # fill with values
    b[:] = T[:] + Q * dt / (rho * c)

    # upper bnd
    if fixed_upper_temperature is not None:
        b[0] = fixed_upper_temperature
    elif upper_bnd_flux is not None:
        b[0] = T[0] + (Q[0] + upper_bnd_flux / dz_top) * dt / (rho[0] * c[0])

    # lower bnd
    if fixed_lower_temperature is not None:
        b[-1] = fixed_lower_temperature
    elif lower_bnd_flux is not None:
        b[-1] = T[-1] + (Q[-1] + lower_bnd_flux / dz_bottom) * dt / (rho[-1] * c[-1])

    return b


def solve_1D_heat_flow_simple(T, dz, dt, K, rho, c, Q,
                              upper_bnd_flux,
                              lower_bnd_flux,
                              fixed_upper_temperature,
                              fixed_lower_temperature,
                              A=None,
                              verbose=False):
    
    """
    solve simple 1D heat flow equation, with constant K, rho, c and dx
    
    specified flux lower bnd condition and specified T upper bnd
    
    """
    if A is None:
        A = construct_heat_flow_matrix(T, dz, dt, K, rho, c,
                                       fixed_upper_temperature,
                                       fixed_lower_temperature)

    nz = len(T)

    b = create_heat_flow_vector(nz, T, Q, dt, rho, c, dz,
                                upper_bnd_flux,
                                lower_bnd_flux,
                                fixed_upper_temperature,
                                fixed_lower_temperature)

    # use numpy linalg to solve system of eqs.:
    # uses  LAPACK routine _gesv:
    # https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/\
    # GUID-4CC4330E-D9A9-4923-A2CE-F60D69EDCAA8.htm
    # uses LU decomposition and iterative solver
    T_new = np.linalg.solve(A, b)

    # TODO check other linear solvers, such as CG:
    # http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/
    # scipy.sparse.linalg.cg.html
    # and use sparse matrices to conserve memory
    # (not really necessary in 1D case)
    
    # check that solution is correct:
    check = np.allclose(np.dot(A, T_new), b)

    if verbose is True:
        print 'solution is correct = ', check

    if check is False:
        print 'warning, solution is correct = ', check

    return T_new, A


def solve_1D_heat_flow(T, z, dt, K, rho, c, Q,
                       upper_bnd_flux,
                       lower_bnd_flux,
                       fixed_upper_temperature,
                       fixed_lower_temperature,
                       A=None,
                       verbose=False):

    """
    solve simple 1D heat flow equation, with variable K, rho, c and dx

    specified flux lower bnd condition and specified T upper bnd

    """

    if A is None:
        A = construct_heat_flow_matrix_variable_z(
            T, z, dt, K, rho, c,
            fixed_upper_temperature,
            fixed_lower_temperature)

    nz = len(T)

    dz_top = z[1] - z[0]
    dz_bottom = z[-1] - z[-2]

    b = create_heat_flow_vector_variable_z(
        nz, T, Q, dt, rho, c, dz_top, dz_bottom,
        upper_bnd_flux,
        lower_bnd_flux,
        fixed_upper_temperature,
        fixed_lower_temperature)

    # use numpy linalg to solve system of eqs.:
    # uses  LAPACK routine _gesv:
    # https://software.intel.com/sites/products/documentation/hpc/mkl/mklman/\
    # GUID-4CC4330E-D9A9-4923-A2CE-F60D69EDCAA8.htm
    # uses LU decomposition and iterative solver
    try:
        T_new = np.linalg.solve(A, b)
    except:
        pdb.set_trace()

    # TODO check other linear solvers, such as CG:
    # http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/
    # scipy.sparse.linalg.cg.html
    # and use sparse matrices to conserve memory
    # (not really necessary in 1D case)

    # check that solution is correct:
    check = np.allclose(np.dot(A, T_new), b)

    if verbose is True:
        print 'solution is correct = ', check

    if check is False:
        print 'warning, solution is ', check

    return T_new, A


def get_geo_history(well_strat, strat_info_mod,
                    max_decompaction_error,
                    exhumation_phase_ids,
                    exhumation_period_starts,
                    exhumation_period_ends,
                    exhumed_strat_units,
                    exhumed_thicknesses,
                    original_thicknesses,
                    max_thickness):

    """
    set up a geohistory dataframe from an input file containing well stratigraphy
    and a file containing compaction parameters of strat units

    :param well_strat:
    :return:
    """

    ############################
    # reconstruct burial history
    ############################
    # deposition code, 1 = burial, -1 = exhumation, 0 = no changes
    well_strat['deposition_code'] = 1

    # get present day thickness
    well_strat['present-day_thickness'] = \
        well_strat['depth_bottom'] - well_strat['depth_top']

    # generate new columns to copy data from strat dataframe
    for col in strat_info_mod.columns:
        well_strat[col] = np.nan

    # copy strat info
    well_strat = copy_df_columns(well_strat, strat_info_mod)

    well_strat = subdivide_strat_units(well_strat, max_thickness)

    # calculate matrix thickness of units (ie thickness with all pore space removed)
    well_strat['matrix_thickness'] = \
        calculate_matrix_thickness(well_strat['surface_porosity'],
                                   well_strat['compressibility'],
                                   well_strat['depth_top'],
                                   well_strat['depth_bottom'])

    # calculate decompacted initial thickness
    print 'decompacting'
    well_strat['initial_thickness'] =\
        calculate_initial_thickness(
            max_decompaction_error,
            well_strat['surface_porosity'],
            well_strat['compressibility'],
            well_strat['depth_top'],
            well_strat['depth_bottom'])

    # add exhumation phases
    geohist_df = \
        add_exhumation_phases(
            well_strat,
            exhumation_phase_ids,
            exhumation_period_starts,
            exhumation_period_ends,
            exhumed_strat_units,
            exhumed_thicknesses,
            original_thicknesses,
            max_thickness,
            strat_info_mod)

    exhumed_units = [g for g in geohist_df.index if g[0] == '+']

    # find age of newly added strat units
    geohist_df = copy_df_columns(geohist_df, strat_info_mod, rows=exhumed_units,
                                 ignore_age_columns=True)

    # and sort again to keep correct order of time
    geohist_df = geohist_df.sort(['age_bottom'])

    # set present-day thickness of fully eroded units to 0
    for strat in geohist_df.index:
        if strat[0] == '+':
            geohist_df.ix[strat, 'present-day_thickness'] = 0
            geohist_df.ix[strat, 'initial_thickness'] = 0

    # calculate matrix thickness of fully eroded units
    for strat in geohist_df.index:
        if strat[0] == '+':
            geohist_df.ix[strat, 'matrix_thickness'] =\
                calculate_matrix_thickness(
                    geohist_df.ix[strat, 'surface_porosity'],
                    geohist_df.ix[strat, 'compressibility'],
                    0,
                    geohist_df.ix[strat, 'eroded_thickness'])

    # find exhumation phases
    exhumation_phases = list(geohist_df.index[geohist_df['deposition_code'] == -1])

    # find maximum burial depth of units:
    geohist_df = find_maximum_depth(geohist_df, exhumation_phases,
                                    max_decompaction_error)

    # remove hiatusses and later eroded units:
    ind = [(s[0] != '~') & (s[0] != '-') for s in geohist_df.index]
    # recalculate matrix thickness using maximum depth instead of
    # present-day depth
    geohist_df['matrix_thickness'][ind] = \
        calculate_matrix_thickness(
            geohist_df['surface_porosity'][ind].values,
            geohist_df['compressibility'][ind].values,
            geohist_df['maximum_depth_top'][ind].values,
            geohist_df['maximum_depth_bottom'][ind].values)

    # find hiatusses
    hiatus_list, hiatus_start_list, hiatus_end_list = \
        find_hiatus(geohist_df.index[::-1],
                    geohist_df['age_bottom'][::-1],
                    geohist_df['age_top'][::-1])

    # store hiatusses in dataframe:
    hiatus_df = pd.DataFrame(index=hiatus_list,
                             columns=geohist_df.columns)
    hiatus_df['deposition_code'] = 0
    hiatus_df['age_bottom'] = hiatus_start_list
    hiatus_df['age_top'] = hiatus_end_list

    # merge hiatusses with geohistory dataframe
    geohist_df = pd.concat((geohist_df, hiatus_df))

    # and sort again to keep correct order of time
    geohist_df = geohist_df.sort(['age_bottom'])

    # set present-day thickness of hiatusses to 0
    for strat in geohist_df.index:
        if strat[0] == '~':
            geohist_df.ix[strat, 'present-day_thickness'] = 0
            geohist_df.ix[strat, 'matrix_thickness'] = 0

    return geohist_df


def reconstruct_strat_thickness(geohist_df):

    """
    create dataframe with thicknesses strat units over time

    :param geohist_df:
    :return:
    """

    # create array of depth vs time:
    times = geohist_df.age_top

    # remove hiatusses and later eroded units:
    ind = [(s[0] != '~') & (s[0] != '-') for s in geohist_df.index]
    strat_units = geohist_df.index[ind]

    strat_thickness_df = pd.DataFrame(index=strat_units,
                                      columns=geohist_df.age_top.astype(str)[::-1])

    # trial forward modeling of compaction
    strat_column = []
    strat_column_list = []
    thickness_list = []

    for timestep in geohist_df.index[::-1]:

        if geohist_df.ix[timestep, 'deposition_code'] == 1:

            # add new strat unit on top
            strat_column.insert(0, timestep)

            # construct lists with compaction params
            n0s = [geohist_df.ix[s, 'surface_porosity'] for s in strat_column]
            betas = [geohist_df.ix[s, 'compressibility'] for s in strat_column]
            bms = [geohist_df.ix[s, 'matrix_thickness'] for s in strat_column]

            # calculate new compacted thicknesses
            thicknesses = np.zeros((len(strat_column)))
            for i, n0, beta, bm in zip(itertools.count(), n0s, betas, bms):
                thicknesses[i] = compact(bm, n0, beta,
                                                    thicknesses[:i].sum(), bm, 0.001)

            if timestep is not geohist_df.index[-1]:
                # check if thickness higher than pre-compacted thickness
                for i, th_old, th_new in zip(itertools.count(), thickness_list[-1], thicknesses[1:]):
                    if th_old < th_new:
                        thicknesses[i+1] = th_old
                #pdb.set_trace()

        # remove top unit
        elif geohist_df.ix[timestep, 'deposition_code'] == -1:
            thicknesses = thickness_list[-1][1:]
            strat_column = strat_column_list[-1][1:]

        #else:
        #    print 'warning, cannot calculate thicknesses'
        #    print 'probably something wrong with the unit ages in the strat_info.csv file'
        #    pdb.set_trace()

        # store strat column and thicknesses in list
        strat_column_list.append(list(strat_column))
        thickness_list.append(thicknesses.copy())

    # now go through burial list and fill the burial history dataframe
    for time, strat_column, thicknesses in zip(times[::-1],
                                               strat_column_list, thickness_list):

        for s, th in zip(strat_column, thicknesses):
            strat_thickness_df.ix[s, str(time)] = th

    # check if modeled and calculated present-day thickness match
    print '\nstrat_unit, thickness, thickness simulated'
    geohist_df['present-day_thickness_simulated'] = np.nan

    lastcol = strat_thickness_df.columns[-1]

    for s, th in zip(strat_column, thicknesses):
        geohist_df.ix[s, 'present-day_thickness_simulated'] = \
            strat_thickness_df.ix[s, lastcol]
        print s, geohist_df.ix[s, 'present-day_thickness'], \
            geohist_df.ix[s, 'present-day_thickness_simulated']

    # fill nan values in thickness df with zeros
    strat_thickness_df = strat_thickness_df.fillna(0.0)

    return strat_thickness_df


def interpolate_param(param_start, param_end, nt_heatflow, fill_empty_value=True):

    """

    :param param_start:
    :param param_end:
    :param nt_heatflow:
    :param fill_empty_value:
    :return:
    """

    int_array = np.resize(np.linspace(0, 1, nt_heatflow +1)[1:],
                          [len(param_start), nt_heatflow]).T

    if fill_empty_value is True:
        a = np.isnan(param_start) == True
        b = np.isnan(param_end) == False
        ind = np.where(a*b == True)[0]
        param_start[ind] = param_end[ind]

        a = np.isnan(param_start) == False
        b = np.isnan(param_end) == True
        ind = np.where(a*b == True)[0]
        param_end[ind] = param_start[ind]

    param_start_array = np.resize(param_start, [nt_heatflow, len(param_start)])
    param_end_array = np.resize(param_end, [nt_heatflow, len(param_start)])

    param_int = int_array * (param_end_array - param_start_array) + param_start_array

    return param_int


def interpolate_node_variables(formation_param_start, formation_param_end,
                               n_timesteps):
    """
    interpolate formation variables to node positions
    node positions: one node at top fm, one at mid and one on bottom

    :param formation_param_start:
    :param formation_param_end:
    :param active_fm_i:
    :param n_timesteps:
    :return:
    """

    formation_var = interpolate_param(formation_param_start,
                                      formation_param_end,
                                      n_timesteps)

    node_var_mid = formation_var.copy()
    node_var_bottom = formation_var.copy()
    node_var_bottom[:, :-1] = (formation_var[:, 1:] +
                               formation_var[:, :-1]) / 2.0
    node_var_bottom[:, -1] = formation_var[:, -1]

    node_var_surface = formation_var[:, 0]

    n_nodes_i = len(formation_param_start) * 2 + 1

    node_var = np.zeros((n_timesteps, n_nodes_i))

    node_var[:, 1::2] = \
        node_var_mid
    node_var[:, 2::2] = \
        node_var_bottom
    node_var[:, 0] = node_var_surface

    return node_var


def generate_thermal_histories(resample_t, n_nodes,
                               time_array_bp,
                               T_nodes, active_nodes,
                               prov_ages_start_array,
                               prov_ages_end_array,
                               nt_prov, Ts):

    aft_node_times = []
    aft_node_temps = []
    #resample_t = pybasin_params.resample_AFT_timesteps

    for nn in xrange(n_nodes):

        # burial history
        burial_time = (time_array_bp[active_nodes[:, nn]][::resample_t] / 1e6)
        burial_T = T_nodes[active_nodes[:, nn], nn][::resample_t]

        T_scenarios = []
        t_scenarios = []

        prov_ages_start_n = prov_ages_start_array[nn]
        prov_ages_end_n = prov_ages_end_array[nn]

        for p1, p2 in zip(prov_ages_start_n,
                          prov_ages_end_n):

            if p1 <= burial_time[0]:
                print 'warning, start deposition in basin earlier than ' \
                      'start provenance history'
                print 'using a hard coded prov history of 120 C to surface'
                print 'from 2 my before deposition to deposition age'
                p1 = burial_time[0] + 2.0
                p2 = burial_time[0]

            if p2 <= burial_time[0]:
                provenance_times = [p1, burial_time[0]]
            else:
                provenance_times = [p1, p2, burial_time[0]]

            # subdivide provenance history into finer timesteps
            prov_time_fine = [np.linspace(a, b, (nt_prov) / (len(provenance_times) - 1) + 1)[:-1]
                              for a, b
                              in zip(provenance_times[:-1],
                                     provenance_times[1:])]

            if len(prov_time_fine) <= 1:
                prov_T_cooling = np.linspace(120.0,
                                             burial_T[0],
                                             len(prov_time_fine[0]))
                prov_T = prov_T_cooling
            else:
                prov_T_depo = np.interp(prov_time_fine[1],
                                        Ts['age'].values,
                                        Ts['surface_temperature'])
                prov_T_cooling = np.linspace(120.0,
                                             prov_T_depo[0],
                                             len(prov_time_fine[0]))
                prov_T = np.concatenate((prov_T_cooling, prov_T_depo))

            prov_time_fine2 = np.concatenate(prov_time_fine)

            combined_times = np.concatenate((prov_time_fine2, burial_time))
            combined_T = np.concatenate((prov_T, burial_T))

            # check if combined time is ok
            if np.max(np.diff(combined_times)) >= 0:
                print 'warning, backwards time in provenance history!'
                pdb.set_trace()

            t_scenarios.append(combined_times.max() - combined_times)
            T_scenarios.append(combined_T)

        aft_node_times.append(t_scenarios)
        aft_node_temps.append(T_scenarios)

    return aft_node_times, aft_node_temps


def generate_burial_histories(n_nodes,
                              time_array_bp,
                              z_nodes, active_nodes,
                              prov_ages_start_array, prov_ages_end_array,
                              nt_prov):

    aft_node_times = []
    aft_node_zs = []
    #resample_t = pybasin_params.resample_AFT_timesteps

    for nn in xrange(n_nodes):

        # burial history
        burial_time = (time_array_bp[active_nodes[:, nn]] / 1e6)
        burial_z = z_nodes[active_nodes[:, nn], nn]

        z_scenarios = []
        t_scenarios = []

        prov_ages_start_n = prov_ages_start_array[nn]
        prov_ages_end_n = prov_ages_end_array[nn]

        for p1, p2 in zip(prov_ages_start_n,
                          prov_ages_end_n):

            if p1 <= burial_time[0]:
                print 'warning, start deposition in basin earlier than ' \
                      'start provenance history'
                print 'using a hard coded prov history of 120 C to surface'
                print 'from 2 my before deposition to deposition age'
                p1 = burial_time[0] + 2.0
                p2 = burial_time[0]

            if p2 <= burial_time[0]:
                provenance_times = [p1, burial_time[0]]
            else:
                provenance_times = [p1, p2, burial_time[0]]

            # subdivide provenance history into finer timesteps
            prov_time_fine = [np.linspace(a, b, (nt_prov) / (len(provenance_times) - 1) + 1)[:-1]
                              for a, b
                              in zip(provenance_times[:-1],
                                     provenance_times[1:])]

            if len(prov_time_fine) <= 1:
                prov_z_cooling = np.linspace(5000.0,
                                             0.0,
                                             len(prov_time_fine[0]))
                prov_z = prov_z_cooling
            else:
                prov_z_depo = np.zeros(len(prov_time_fine[1]))
                prov_z_cooling = np.linspace(5000.0,
                                             0.0,
                                             len(prov_time_fine[0]))
                prov_z = np.concatenate((prov_z_cooling, prov_z_depo))

            prov_time_fine2 = np.concatenate(prov_time_fine)

            combined_times = np.concatenate((prov_time_fine2, burial_time))
            combined_z = np.concatenate((prov_z, burial_z))

            # check if combined time is ok
            if np.max(np.diff(combined_times)) >= 0:
                print 'warning, backwards time in provenance history!'
                pdb.set_trace()

            t_scenarios.append(combined_times)
            z_scenarios.append(combined_z)

        aft_node_times.append(t_scenarios)
        aft_node_zs.append(z_scenarios)

    return aft_node_times, aft_node_zs


def calculate_aft_ages_pdf(aft_ages, aft_ages_min_std, aft_ages_plus_std):
    """
    calculate pdf of fission track ages
    convert age to gamma (Galbraith ....), assume normal distribution of gamma
    and use this to calculate aft age distribution pdf

    :param aft_ages:
    :param aft_ages_min_std:
    :param aft_ages_plus_std:
    :return:
    """

    binrange = 30
    binsize = 0.01
    bins = np.arange(0, binrange, binsize)
    bins[0] = 1.0e-10

    aft_ages_min = aft_ages - aft_ages_min_std
    aft_ages_max = aft_ages + aft_ages_plus_std

    # calculate gamma value (i.e. normally distributed transform of AFT age equation)
    gamma = AFTlib.calculate_gamma_from_AFT_age(aft_ages.values)
    #print 'gamma :', gamma

    # correct zero ages:
    zero_index = np.where(aft_ages == 0)[0]
    gamma[zero_index] = 1e-10

    # calculate standard error of gamma:
    gamma_max = AFTlib.calculate_gamma_from_AFT_age(aft_ages_max.values)
    gamma_std = (gamma_max - gamma)/2.0

    gamma_pdf = np.zeros((len(aft_ages), len(bins)))

    for i in xrange(len(aft_ages)):
        gamma_pdf[i] = scipy.stats.norm.pdf(bins, gamma[i], gamma_std[i])

        # normalize pdf to a value of 1:
        gamma_pdf[i] = gamma_pdf[i]/(gamma_pdf[i].sum())

    #return gamma_pdf
    aft_age_bins = AFTlib.calculate_AFT_age_from_gamma(bins)

    return aft_age_bins, gamma_pdf


def calculate_vr(T_nodes, active_nodes, time_array, n_nodes):

    vr_nodes = np.zeros(T_nodes.shape)
    for nn in xrange(n_nodes):
        vr_nodes[active_nodes[:, nn], nn] = \
            easyRo.easyRo(time_array[active_nodes[:, nn]] / 1e6,
                          T_nodes[active_nodes[:, nn], nn])

    return vr_nodes


def simulate_aft(resample_t, nt_prov, n_nodes, time_array_bp, z_nodes, T_nodes, active_nodes,
                 prov_ages_start, prov_ages_end, annealing_kinetics_values, annealing_kinetic_param, Ts):

    """
    simulate fission track ages using calculated burial thermal history and provenance thermal history scenarios

    :param resample_t:
    :param nt_prov:
    :param n_nodes:
    :param time_array_bp:
    :param z_nodes:
    :param T_nodes:
    :param active_nodes:
    :param prov_ages_start:
    :param prov_ages_end:
    :param annealing_kinetics_values:
    :param annealing_kinetic_param:
    :param Ts:
    :return:
    """

    # combine burial and provenance history

    aft_node_times, aft_node_temps = generate_thermal_histories(
        resample_t, n_nodes,
        time_array_bp, T_nodes, active_nodes,
        prov_ages_start, prov_ages_end,
        nt_prov, Ts)

    aft_node_times_burial, aft_node_zs = generate_burial_histories(
        n_nodes,
        time_array_bp, z_nodes, active_nodes,
        prov_ages_start, prov_ages_end,
        nt_prov)

    #figb = pl.figure()
    #ax = figb.add_subplot(1, 1, 1)
    #for x, y in zip(aft_node_times_burial, aft_node_zs):
    #    for xi, yi in zip(x, y):
    #        ax.plot(xi, yi)
    #figb.savefig('burial_prov_test.png', dpi=300)

    # calculate FT ages for all formations
    print 'calculating FT ages and lengths for all nodes'
    n_prov_scenarios = prov_ages_start.shape[1]
    n_kinetic_scenarios = len(annealing_kinetics_values)
    aft_age_nodes = np.zeros((n_nodes, n_prov_scenarios,
                              n_kinetic_scenarios))
    aft_ln_mean_nodes = np.zeros((n_nodes, n_prov_scenarios,
                                  n_kinetic_scenarios))
    aft_ln_std_nodes = np.zeros((n_nodes, n_prov_scenarios,
                                 n_kinetic_scenarios))

    for nn in xrange(n_nodes):

        for n_prov in xrange(n_prov_scenarios):
            for n_kin in xrange(n_kinetic_scenarios):
                trackLengthPDF, AFTage, l_mean, l_mean_std, rm, rc, rho_age, dt =\
                    AFTannealingLib.simulate_AFT_annealing(
                        aft_node_times[nn][n_prov],
                        aft_node_temps[nn][n_prov],
                        annealing_kinetics_values[n_kin],
                        kinetic_parameter=annealing_kinetic_param,
                        use_fortran_algorithm=True)

                aft_age_nodes[nn, n_prov, n_kin] = AFTage
                aft_ln_mean_nodes[nn, n_prov, n_kin] = l_mean
                aft_ln_std_nodes[nn, n_prov, n_kin] = l_mean_std

    aft_age_nodes_min = np.min(aft_age_nodes, axis=(1, 2))
    aft_age_nodes_max = np.max(aft_age_nodes, axis=(1, 2))

    return (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
            aft_ln_mean_nodes, aft_ln_std_nodes,
            aft_node_times_burial, aft_node_zs)


def run_burial_hist_model(well_number, well, well_strat, strat_info_mod,
                          pybasin_params,
                          Ts, litho_props,
                          output_dir, model_scenario_number):
    """
    run burial and thermal history model

    :param well_number:
    :param well:
    :param well_strat:
    :param strat_info_mod:
    :param pybasin_params:
    :param output_dir:
    :param Ts:
    :param litho_props:
    :return:
    """
    year = 365.25 * 24.0 * 60 * 60.

    geohist_df = get_geo_history(
        well_strat, strat_info_mod,
        pybasin_params.max_decompaction_error,
        pybasin_params.exhumation_phase_ids,
        pybasin_params.exhumation_period_starts,
        pybasin_params.exhumation_period_ends,
        pybasin_params.exhumed_strat_units,
        pybasin_params.exhumed_thicknesses,
        pybasin_params.original_thicknesses,
        pybasin_params.max_thickness)

    # save geohistory dataframe as .csv file
    geohist_df.to_csv(os.path.join(output_dir,
                                   'geo_history_well_%s_ms%i.csv'
                                   % (well, model_scenario_number)),
                      index_label='strat_unit')

    strat_thickness_df = \
        reconstruct_strat_thickness(geohist_df)

    strat_thickness_df.to_csv(
        os.path.join(output_dir,
                     'formation_thickness_history_well_%s_ms%i.csv'
                     % (well, model_scenario_number)),
        index_label='strat_unit')

    # calculate burial depths from thickness dataframe
    burial_df = strat_thickness_df.copy()
    for ind in burial_df.index:
        burial_df.ix[ind] = strat_thickness_df.ix[:ind].sum()

    burial_df.to_csv(os.path.join(output_dir,
                                  'burial_history_well_%s_ms%i.csv'
                                  % (well, model_scenario_number)),
                     index_label='strat_unit')

    # find which formations are active and inactive
    active_fm = pd.DataFrame(index=strat_thickness_df.index,
                             columns=strat_thickness_df.columns,
                             dtype=bool)

    for start_age, end_age in zip(strat_thickness_df.columns[:-1],
                                  strat_thickness_df.columns[1:]):
        inactive_bool = np.logical_and(strat_thickness_df[start_age]==0,
                                       strat_thickness_df[end_age]==0)
        active_bool = inactive_bool.values == False

        try:
            active_fm_backup = active_fm.copy()
            active_fm[start_age] = active_bool
        except TypeError, msg:
            print msg
            pdb.set_trace()

    active_fm[end_age] = active_fm[start_age]

    # get ages
    ages_ind = burial_df.columns
    ages = burial_df.columns.values.astype(float)

    # calculate duration of geological timesteps
    durations = -np.diff(ages)

    # populate arrays with thermal conductivity, heat capacity, heat production,
    # density and matrix thickness
    hf_param_df = pd.DataFrame(index=burial_df.index,
                               columns=['K', 'c', 'HP',
                                        'rho_matrix',
                                        'matrix_thickness',
                                        'surface_porosity'])

    for ix in hf_param_df.index:
        hf_param_df.ix[ix, 'K'] = geohist_df.ix[ix, 'thermal_conductivity']
        hf_param_df.ix[ix, 'c'] = geohist_df.ix[ix, 'heat_capacity']
        hf_param_df.ix[ix, 'HP'] = geohist_df.ix[ix, 'heat_production']
        hf_param_df.ix[ix, 'rho_matrix'] = 2650.0
        hf_param_df.ix[ix, 'matrix_thickness'] = geohist_df.ix[ix, 'matrix_thickness']
        hf_param_df.ix[ix, 'surface_porosity'] = geohist_df.ix[ix, 'surface_porosity']

    # calculate evolution of porosity
    porosity_df = pd.DataFrame(index=burial_df.index,
                               columns=burial_df.columns)
    k_df = pd.DataFrame(index=burial_df.index,
                        columns=burial_df.columns)
    c_df = pd.DataFrame(index=burial_df.index,
                        columns=burial_df.columns)
    rho_df = pd.DataFrame(index=burial_df.index,
                          columns=burial_df.columns)
    hp_df = pd.DataFrame(index=burial_df.index,
                         columns=burial_df.columns)

    # calculate porosity at each timestep
    for ix in porosity_df.index:
        porosity_df.ix[ix] = (1.0
                              - hf_param_df.ix[ix, 'matrix_thickness']
                              / strat_thickness_df.ix[ix])

    # change inf and -inf to nan values
    porosity_df[porosity_df == -np.inf] = np.nan
    porosity_df[porosity_df == np.inf] = np.nan

    # set porosity of active formations with 0 initial thickness to
    # surface porosity
    for col in porosity_df.columns:
        for row in porosity_df.index:
            if (active_fm.ix[row, col] == True
                    and np.isnan(porosity_df.ix[row, col])):
                porosity_df.ix[row, col] = hf_param_df.ix[row,
                                                          'surface_porosity']

    # calculate bulk thermal conductivity, heat capacity, density and
    # heat production
    # at each geological timestep
    for ix in c_df.index:
        k_df.ix[ix] = (hf_param_df.ix[ix, 'K'] ** (1.0 - porosity_df.ix[ix]) *
                       litho_props.ix['water', 'thermal_conductivity'] ** porosity_df.ix[ix])
        c_df.ix[ix] = (hf_param_df.ix[ix, 'c'] * (1.0 - porosity_df.ix[ix]) +
                       litho_props.ix['water', 'heat_capacity'] * porosity_df.ix[ix])
        rho_df.ix[ix] = (hf_param_df.ix[ix, 'rho_matrix'] *
                         (1.0 - porosity_df.ix[ix]) +
                         porosity_df.ix[ix] *
                         litho_props.ix['water', 'density'])
        hp_df.ix[ix] = hf_param_df.ix[ix, 'HP'] * (1.0 - porosity_df.ix[ix])

    ##############################################
    # set up arrays for forward model of heat flow
    ##############################################

    # calculate timesteps
    nt_heatflows = np.array([int(np.round(duration * 1e6 /
                                          pybasin_params.max_hf_timestep))
                             for duration in durations])

    nt_total = nt_heatflows.sum()

    dt_hfs = durations * 1.0e6 / nt_heatflows

    # grid cell points at top, mid and bottom of formations
    n_cells = burial_df.shape[0] * 2
    n_nodes = n_cells + 1

    # set up arrays for depth, temperature and thermal parameters
    z_nodes = np.zeros((nt_total, n_nodes))
    node_age = np.zeros(n_nodes)
    k_nodes = np.zeros((nt_total, n_cells))
    c_nodes = np.zeros((nt_total, n_nodes))
    rho_nodes = np.zeros((nt_total, n_nodes))
    hp_nodes = np.zeros((nt_total, n_nodes))

    start_ages = ages_ind[:-1]
    end_ages = ages_ind[1:]

    # find top formation at each timestep:
    top_fms = [np.where(active_fm[start_age] == True)[0][0]
               for start_age in start_ages]

    # find out which formations were active
    active_cells = np.zeros((nt_total, n_cells), dtype=bool)
    active_nodes = np.zeros((nt_total, n_cells + 1), dtype=bool)
    timestep = 0
    for start_age, end_age, nt_heatflow in \
            zip(start_ages, end_ages, nt_heatflows):

        active_fm_i = active_fm[start_age].values

        active_cells[timestep:(timestep + nt_heatflow), :-1:2] = \
            active_fm_i
        active_cells[timestep:(timestep + nt_heatflow), 1::2] = \
            active_fm_i

        active_nodes[timestep:(timestep + nt_heatflow), 1:-1:2] = \
            active_fm_i
        active_nodes[timestep:(timestep + nt_heatflow), 2::2] = \
            active_fm_i

        # add surface node
        i = np.where(active_fm_i == True)[0][0]
        active_nodes[timestep:(timestep + nt_heatflow), i*2] = True

        timestep += nt_heatflow

    # calculate depths:
    timestep = 0
    for start_age, end_age, nt_heatflow in \
            zip(start_ages, end_ages, nt_heatflows):


        try:
            zs_bottom = interpolate_param(burial_df[start_age].values,
                                          burial_df[end_age].values,
                                          nt_heatflow)
        except ValueError, msg:
            print msg
            pdb.set_trace()
        zs_surface = np.zeros((nt_heatflow, 1))
        zs_top = np.hstack((zs_surface, zs_bottom[:, :-1]))
        #zs_all = np.hstack((zs_surface, zs_bottom))

        zs_mid = (zs_bottom + zs_top) / 2.0

        # add thicknesses to final array
        z_nodes[timestep:(timestep + nt_heatflow), 1::2] = zs_mid
        z_nodes[timestep:(timestep + nt_heatflow), 2::2] = zs_bottom

        timestep += nt_heatflow

    thickness_all = np.diff(z_nodes, axis=1)

    node_strat_r = burial_df.index.tolist()
    node_strat = [''] * n_nodes
    node_strat[1::2] = node_strat_r
    node_strat[2::2] = node_strat_r
    node_strat[0] = node_strat[1]

    cell_strat_r = burial_df.index.tolist()
    cell_strat = [''] * n_cells
    cell_strat[::2] = cell_strat_r
    cell_strat[1::2] = cell_strat_r
    cell_strat[0] = cell_strat[1]

    ind = geohist_df['deposition_code'] == 1
    age_nodes_start_raw = geohist_df['age_bottom'][ind]
    age_nodes_end_raw = geohist_df['age_top'][ind]

    node_age[:-1:2] = age_nodes_end_raw
    node_age[-1] = age_nodes_start_raw[-1]
    node_age[1::2] = (age_nodes_end_raw + age_nodes_start_raw) / 2.0

    # copy provnenacne ages for nodes
    prov_cols = [col for col in geohist_df.columns if 'provenance' in col]
    n_prov = len(prov_cols) / 2
    prov_ages_start = []
    prov_ages_end = []

    for i in xrange(n_prov):
        prov_ages_start.append(geohist_df['provenance_age_start_%i'
                                          % (i + 1)])
        prov_ages_end.append(geohist_df['provenance_age_end_%i'
                                        % (i + 1)])

    prov_ages_start_array = np.array(prov_ages_start).T
    prov_ages_end_array = np.array(prov_ages_end).T

    # filter out hiatusses
    active_prov = np.isnan(prov_ages_start_array[:, -1].astype(float)) == False

    prov_start_nodes = np.zeros((n_nodes, n_prov))
    prov_end_nodes = np.zeros((n_nodes, n_prov))

    prov_start_nodes[1::2] = prov_ages_start_array[active_prov]
    prov_start_nodes[2::2] = prov_ages_start_array[active_prov]
    prov_start_nodes[0] = prov_start_nodes[1]
    
    prov_end_nodes[1::2] = prov_ages_end_array[active_prov]
    prov_end_nodes[2::2] = prov_ages_end_array[active_prov]
    prov_end_nodes[0] = prov_end_nodes[1]

    # populate thermal parameter arrays:
    print 'setting grid cell thermal params (K)'
    timestep = 0
    for start_age, end_age, nt_heatflow in \
            zip(start_ages, end_ages, nt_heatflows):

        ks = interpolate_param(k_df[start_age].values.astype(float),
                               k_df[end_age].values.astype(float),
                               nt_heatflow)

        k_nodes[timestep:(timestep + nt_heatflow), :-1:2] = ks
        k_nodes[timestep:(timestep + nt_heatflow), 1::2] = ks

        timestep += nt_heatflow

    # populate thermal parameter arrays:
    print 'setting grid node thermal params (rho, c, HP)'
    timestep = 0
    for start_age, end_age, nt_heatflow in \
            zip(start_ages, end_ages, nt_heatflows):

        end_step = (timestep + nt_heatflow)

        active_fm_i = active_fm[start_age].values

        rho_active = interpolate_node_variables(
            rho_df[start_age].values.astype(float)[active_fm_i],
            rho_df[end_age].values.astype(float)[active_fm_i],
            nt_heatflow)

        rho_nodes[timestep:end_step, active_nodes[timestep]] = rho_active

        c_active = interpolate_node_variables(
            c_df[start_age].values.astype(float)[active_fm_i],
            c_df[end_age].values.astype(float)[active_fm_i],
            nt_heatflow)

        c_nodes[timestep:end_step, active_nodes[timestep]] = c_active

        hp_active = interpolate_node_variables(
            hp_df[start_age].values.astype(float)[active_fm_i],
            hp_df[end_age].values.astype(float)[active_fm_i],
            nt_heatflow)

        hp_nodes[timestep:end_step, active_nodes[timestep]] = hp_active

        timestep += nt_heatflow

    # check to remove 0 depth nodes after erosion phases
    for nti in xrange(nt_total):

        ind0 = np.where(np.diff(
            z_nodes[nti, active_nodes[nti]])==0)[0]

        if len(ind0) > 0:
            active_nodes[nti,
                         np.diff(z_nodes[nti])==0] = False
            active_cells[nti] = active_nodes[nti, :-1]

    # check grid nodes
    active_grid_nodes_sum = np.sum(active_nodes, axis=1)

    # generate time array
    time_array = np.zeros(nt_total)
    timestep = 0
    time_all = 0
    for dt_hf, nt_heatflow in zip(dt_hfs, nt_heatflows):

        time_array[timestep:(timestep + nt_heatflow)] = \
            np.arange(nt_heatflow) * dt_hf + time_all
        time_all += nt_heatflow * dt_hf
        timestep += nt_heatflow

    # calculate time in yr bp
    time_array_bp = time_array.max() - time_array

    # interpolate surface temperature
    surface_temp_array = np.interp(time_array_bp,
                                   Ts['age'].values * 1.0e6,
                                   Ts['surface_temperature'])

    # set initial temp
    T_nodes = np.zeros((nt_total, n_nodes))
    for ni in xrange(nt_total):
        T_nodes[ni, :] = surface_temp_array[ni]

    # interpolate basal heat flow
    basal_hf_array = np.interp(time_array_bp,
                               pybasin_params.heatflow_ages * 1.0e6,
                               pybasin_params.heatflow_history)

    # go through all geological timesteps and model heat flow:
    print '-' * 10
    print 'modeling heatflow'

    cumulative_steps = np.cumsum(nt_heatflows)

    for timestep in xrange(nt_total):

        active_cells_i = active_cells[timestep]
        active_nodes_i = active_nodes[timestep]

        if timestep == 0:
            T_init = T_nodes[timestep, active_nodes_i]
        else:
            T_init = T_nodes[timestep-1, active_nodes_i]

        if np.any(np.isnan(T_init)):
            print 'error'
            print bla

        # calculate temperature
        T_nodes[timestep, active_nodes_i], A = \
            solve_1D_heat_flow(
                T_init,
                z_nodes[timestep, active_nodes_i],
                dt_hf * year,
                k_nodes[timestep, active_cells_i],
                rho_nodes[timestep, active_nodes_i],
                c_nodes[timestep, active_nodes_i],
                hp_nodes[timestep, active_nodes_i],
                None,
                basal_hf_array[timestep],
                surface_temp_array[timestep],
                None)

        if np.any(np.isnan(T_nodes[timestep, active_nodes_i])):
            print 'error, nan values in T array'

            for cs, ac, rho in zip(cell_strat, active_cells[timestep],
                                   rho_nodes[timestep]):

                print cs, ac, rho

            print bla

        if timestep in cumulative_steps:

            age_ind = np.where(cumulative_steps==timestep)[0][0]

            print 'step %i, %0.2f Ma, max z = %0.1f, min, max T = %0.1f - %0.1f' \
                  % (timestep, time_array_bp[timestep] / 1e6,
                     z_nodes[timestep, active_nodes_i].max(),
                     T_nodes[timestep, active_nodes_i].min(),
                     T_nodes[timestep, active_nodes_i].max())

    return (geohist_df, time_array, time_array_bp,
            surface_temp_array, basal_hf_array,
            z_nodes, T_nodes, active_nodes,
            n_nodes, n_cells,
            node_strat, node_age,
            prov_start_nodes, prov_end_nodes)