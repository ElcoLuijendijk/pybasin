"""
functions for reading and checking pybasin model input data

this is the input-side counterpart to lib/pybasin_io.py, which handles
saving and loading model run *output* data. the functions here check
that a model input directory has all the required .csv files, read
those files into pandas DataFrames, and select the subset of input
data (well stratigraphy, surface salinity boundary condition) that
applies to one particular well
"""

import logging
import os
import re

import numpy as np
import pandas as pd

from . import pybasin_lib

logger = logging.getLogger(__name__)


def check_input_data_files(input_dir, pybasin_params):

    """
    check if all necessary input files are available

    """

    logger.info('checking for input files in %s' % input_dir)

    fns = ['stratigraphy_info.csv', 'well_stratigraphy.csv', 'surface_temperature.csv',
           'lithology_properties.csv', 'temperature_data.csv']

    if pybasin_params.simulate_salinity is True:
        fns += ['surface_salinity.csv', 'salinity_data.csv']

    if pybasin_params.simulate_VR is True:
        fns.append('vitrinite_reflectance.csv')

    if pybasin_params.simulate_He is True:
        fns += ['he_samples.csv', 'he_data.csv']

    if pybasin_params.simulate_AFT is True:
        fns += ['aft_samples.csv', 'aft_data.csv']

    for fn in fns:
        if os.path.exists(os.path.join(input_dir, fn)) is False:
            msg = 'could not find input file %s in input directory %s' % (fn, input_dir)
            raise FileNotFoundError(msg)

    logger.info('found all necessary input files in %s' % input_dir)

    return


def read_model_input_data(input_dir, pybasin_params):

    """
    Read all model input data .csv files using Pandas

    """

    # read all input data
    logger.info('reading input data')
    # stratigraphy description
    strat_info = pd.read_csv(os.path.join(input_dir, 'stratigraphy_info.csv'), skip_blank_lines=True)


    required_cols = ['strat_unit', 'age_bottom', 'age_top']

    for required_col in required_cols:
        try:
            assert required_col in strat_info.columns
        except AssertionError:
            fn_msg = os.path.join(input_dir, 'stratigraphy_info.csv')
            msg = f"could not find the column {required_col} in the stratigraphy input file {fn_msg}\n"
            msg += f"instead I found only the following {len(strat_info.columns)} columns: {list(strat_info.columns)}\n"
            msg += f"please check the file with a text editor to see if {required_col} is present in the header"
            msg += f"and if the column separators are consistent"
            raise ValueError(msg)

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

    # for the optional effective stress-based compaction method: make
    # sure a compressibility_stress column (Pa^-1) is present, and
    # auto-derive it from the existing depth-based compaction
    # coefficient (m^-1, stored in the 'compressibility' column) for
    # any lithology that does not have an explicit value. this means
    # existing lithology_properties.csv files keep working unchanged,
    # since this column and the auto-derivation are only used if
    # compaction_method is set to 'effective_stress'
    compaction_method = getattr(pybasin_params, 'compaction_method', 'depth')
    if compaction_method == 'effective_stress':
        if 'compressibility_stress' not in litho_props.columns:
            litho_props['compressibility_stress'] = np.nan

        water_density = litho_props.loc['water', 'density']

        for lith in litho_props.index:
            if lith == 'water':
                continue
            if pd.isnull(litho_props.loc[lith, 'compressibility_stress']):
                density_contrast = litho_props.loc[lith, 'density'] - water_density
                if density_contrast <= 0:
                    msg = (f"cannot auto-derive compressibility_stress for lithology "
                           f"'{lith}': its density ({litho_props.loc[lith, 'density']}) is "
                           f"not greater than the water density ({water_density}). Please "
                           f"provide an explicit compressibility_stress value for this "
                           f"lithology in lithology_properties.csv")
                    raise ValueError(msg)
                litho_props.loc[lith, 'compressibility_stress'] = (
                    litho_props.loc[lith, 'compressibility']
                    / (density_contrast * pybasin_lib.G_ACCEL))

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

    if pybasin_params.simulate_He is True:
        he_sample_fn = os.path.join(input_dir, 'he_samples.csv')
        he_data_fn = os.path.join(input_dir, 'he_data.csv')
        if os.path.exists(he_sample_fn):
            # read U-Th/He (He) data
            he_samples = pd.read_csv(he_sample_fn, skip_blank_lines=True)
        else:
            logger.warning('could not find input file %s' % he_sample_fn)
            logger.info('continuing without U-Th/He sample data')
            he_samples = None

        if os.path.exists(he_data_fn):
            # read apatite U-Th/He (He) data
            he_data = pd.read_csv(he_data_fn, skip_blank_lines=True)
        else:
            logger.warning('could not find input file %s' % he_data_fn)
            logger.info('continuing without  U-Th/He age data')
            he_data = None

    else:
        he_samples = None
        he_data = None

    if pybasin_params.simulate_salinity is True:
        # read surface salinity bnd condditions
        # Cs = pd.read_csv(os.path.join(input_dir, 'surface_salinity.csv'))

        # and read salinity data
        salinity_data = pd.read_csv(os.path.join(input_dir, 'salinity_data.csv'), skip_blank_lines=True)

    else:
        salinity_data = None

    # optional observed present-day pressure data (eg. from drill stem
    # tests), used both to overlay on the pressure panel of the burial
    # history figure and to calculate the goodness of fit of modeled
    # vs measured fluid pressure
    pressure_data_file = os.path.join(input_dir, 'dst_pressure_data.csv')
    if getattr(pybasin_params, 'simulate_fluid_flow', False) is True \
            and os.path.isfile(pressure_data_file):
        pressure_data = pd.read_csv(pressure_data_file, skip_blank_lines=True)
    else:
        pressure_data = None

    # optional observed porosity data, used only to overlay on the
    # porosity panel of the burial history figure. not used for
    # calibration or goodness-of-fit
    porosity_data_file = os.path.join(input_dir, 'porosity_data.csv')
    if os.path.isfile(porosity_data_file):
        porosity_data = pd.read_csv(porosity_data_file, skip_blank_lines=True)
    else:
        porosity_data = None

    # T_data, vr_data, aft_samples, aft_ages, he_samples, he_data, salinity_data

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
    except AssertionError as e:
        msg = ('not all lithology units found in the stratigraphy info file '
              'are also present in lithology_properties.csv: %s' % e)
        logger.error(msg)
        raise AssertionError(msg) from e

    # provenance histories are numbered starting at 1 (provenance_age_start_1,
    # provenance_age_end_1, provenance_age_start_2, ...). the highest-numbered
    # histories are optional: drop any that are entirely unused, as long as
    # at least min_provenance_histories remain with data
    min_provenance_histories = 2

    prov_history_numbers = sorted(set(
        int(re.search(r'_(\d+)$', col).group(1)) for col in prov_cols))

    cols_to_drop = []
    for n in sorted(prov_history_numbers, reverse=True):
        start_col = 'provenance_age_start_%i' % n
        end_col = 'provenance_age_end_%i' % n
        history_is_empty = bool(strat_info[start_col].isnull().all()
                                 and strat_info[end_col].isnull().all())

        if not history_is_empty:
            break

        cols_to_drop += [start_col, end_col]

    n_remaining_histories = len(prov_history_numbers) - len(cols_to_drop) // 2

    if n_remaining_histories < min_provenance_histories:
        msg = ("error in parsing stratigraphy_info.csv: found data for only "
              "%i provenance history(s), at least %i provenance histories "
              "are required" % (n_remaining_histories, min_provenance_histories))
        raise ValueError(msg)

    if cols_to_drop:
        logger.info('stratigraphy_info.csv contains %i unused provenance '
                   'history column(s) (%s), these will be ignored'
                   % (len(cols_to_drop) // 2, ', '.join(cols_to_drop)))
        strat_info = strat_info.drop(columns=cols_to_drop)
        prov_cols = [col for col in prov_cols if col not in cols_to_drop]

    # check remaining provenance columns for unexpected gaps, ie. an empty
    # column followed by a populated higher-numbered one
    for p in prov_cols:
        if np.all(strat_info[p].isnull()):
            msg = ("error in parsing stratigraphy_info.csv: the provenance age "
                  "column '%s' is empty for every stratigraphic unit" % p)
            raise ValueError(msg)

    # check well stratigraphy file
    if well_strats['depth_top'].min() < 0 or well_strats['depth_bottom'].min() < 0:

        msg = ('found a negative value for depth in the depth_top or '
              'depth_bottom columns of well_stratigraphy.csv. Please make sure all values for '
              'depth are zero or positive')

        raise ValueError(msg)

    # create new copy of dataframe to store results
    strat_info_mod = strat_info.copy()

    # lithology properties that are not fraction-weighted blended into
    # strat_info_mod, either because they are non-numeric (eg. a
    # mineral name) or because the model reads them directly from
    # litho_props per pure lithology instead (eg. permeability
    # end-member inputs used by calculate_bulk_permeability_df())
    non_blended_litho_props = ['clay_mineral']

    numeric_litho_props = [col for col in litho_props.columns
                           if col not in non_blended_litho_props]

    # add new lithology properties columns to stratigraphy dataframe
    for litho_prop_name in numeric_litho_props:
        strat_info_mod[litho_prop_name] = 0

    # go through all litho properties
    for litho_prop_name in numeric_litho_props:

        # go through all lithology types present in strat unit
        for col in litho_cols:
            # find fraction of lithology and multiply by lithology property
            s = strat_info[col].astype(float) * litho_props[litho_prop_name][col]

            # add fraction to column in strat dataframe
            strat_info_mod[litho_prop_name] = strat_info_mod[litho_prop_name] + s

    return (well_strats, strat_info_mod, df_sal,
            T_data_df, vr_data_df,
            aft_samples, aft_ages,
            he_samples, he_data,
            salinity_data, surface_temp, litho_props,
            pressure_data, porosity_data)


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
        raise ValueError(msg)

    return well_strat, well_strat_orig


def select_well_salinity_bnd(well, salinity_bnd_df):

    """
    Read groundwater salinity values at top boundary
    """

    if well in salinity_bnd_df.columns:
        logger.info('using surface salinity bnd for well %s from file' % well)
        # cols = ['age_start', 'age_end', 'surface_salinity_%s' % well]
        cols = ['age_start', 'age_end', well]
        surface_salinity_well = salinity_bnd_df[cols]
        surface_salinity_well['surface_salinity'] = \
            surface_salinity_well[well]
    else:
        surface_salinity_well = None

    return surface_salinity_well
