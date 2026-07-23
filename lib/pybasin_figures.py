"""
module that contains all functions for making figures of pybasin model results

"""


__author__ = 'elcopone'

import itertools
import math
import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
#import matplotlib.mlab
import matplotlib.patches as mpatches
from matplotlib import ticker
import matplotlib

from . import pybasin_io

logger = logging.getLogger(__name__)


def setup_figure(width=125.0, height='g', fontsize='x-small',
                 landscape=False, fontsize_legend=0, units='mm'):

    """
    Set up a new Matplotlib figure

    Default figure is 125 mm wide and 125 mm/ golden ratio high

    figure size AGU: one column = 8.4 cm, 2 column = 16.98 cm, max height = 23.7

    Parameters
    ----------
    width : float, optional
        horizontal size of image, default is 125.0 (mm)
        use '1col' and '2col' for default widths agu journal figs
        (84.0 and 169.8 mm)
    height : float, optional
        vertical size of figure
        if 'g' the height is equal to the width divided by the golden ration
        if floating point number, height = width * number
    fontsize : string, optional
        'xxx-small', 'xx-small', 'x-small', 'small' or 'medium'
        default is 'x-small'
    landscape : boolean, optional
        Landscape figure instead of portrait
        default is False
    fontsize_legend : float or string, optional
        text size of the legend, use default matplotlib fontsize formatting
        if fotsize is None the fontsize of legend items is determined by the
        default settings of the fontsize parameter
        default is None
    units : string, optional
        units used to determine figure size, choose either 'inch' or 'mm'
        default = 'mm'

    Returns
    -------
    fig : matplotlib figure instance
    """

    if width == '1col':
        width = 84.0
    elif width == '2col':
        width = 169.8

    golden_ratio = (1.0 + np.sqrt(5))/2.0

    if type(height) == str and height[-1] == 'g':
        if len(height) > 1:
            c = float(height[:-1])
        else:
            c = 1
        height = width / (c * golden_ratio)
    elif type(height) == float or type(height) == int:
        height = width * height
    elif height == 'max':
        height = 170
    if height > 215.0:
        logger.info('figure exceeding b5 paper size')

    # initialize figure
    if landscape is True:
        logger.info('landscape figure')
        xs = height
        ys = width
    else:
        logger.info('portrait figure')
        xs = width
        ys = height

    if units != 'inch':
        xs = xs / 25.4
        ys = ys / 25.4

    logger.info('init fig,  size = %0.1f x %0.1f inch' % (xs, ys))

    fig = pl.figure(figsize=(xs,ys))

    # set default parameters for figure
    if type(fontsize) == str:
        if fontsize == 'xxx-small':
            fontsize_s = 'xx-small'
            fontsize_l = 'xx-small'
            fontsize_leg = 'xx-small'
        elif fontsize == 'xx-small':
            fontsize_s = 'xx-small'
            fontsize_l = 'x-small'
            fontsize_leg = 'xx-small'
        elif fontsize == 'x-small':
            fontsize_s = 'x-small'
            fontsize_l = 'small'
            fontsize_leg = 'xx-small'
        elif fontsize == 'small':
            fontsize_s = 'small'
            fontsize_l = 'medium'
            fontsize_leg = 'x-small'
        elif fontsize == 'medium':
            fontsize_s = 'medium'
            fontsize_l = 'large'
            fontsize_leg = 'small'
    else:
        fontsize_s = fontsize
        fontsize_l = 'xx-small'
        fontsize_leg = fontsize

    if fontsize_legend is not None:
        fontsize_leg = fontsize_legend

    params = {'axes.labelsize': fontsize_s,
              'text.fontsize': fontsize_l,
              'legend.fontsize': fontsize_leg,
              'axes.titlesize': fontsize_l,
              'xtick.labelsize': fontsize_s,
              'ytick.labelsize': fontsize_s}

    #pl.rcParams.update(params)

    return fig


def model_vs_data_figure(model_run_data,
                         show_provenance_hist=True,
                         show_strat_column=False,
                         show_thermochron_data=True,
                         contour_variable='temperature',
                         add_legend=True,
                         strat_fontsize='xx-small',
                         figsize=170.0,
                         legend_space=0.12,
                         height_ratio=3,
                         cb_buffer_vert=-0.03,
                         cb_buffer_hor=0.0,
                         show_prov_ages_simple=False,
                         show_violin_plot=True,
                         ncols_legend=3,
                         add_panel_titles=True,
                         panel_title_numbers=False,
                         panel_title_prefix='',
                         panel_label_fs='small',
                         legend_fontsize='x-small',
                         bottom=0.12,
                         left=0.12,
                         right=0.97,
                         top=0.96,
                         max_strat_units=None,
                         max_age_burial_panel=None,
                         max_age_thermochron_panel=None,
                         show_temp_panel=True,
                         show_porosity_panel=None,
                         show_pressure_panel=None,
                         show_age_scatter=True,
                         show_gof_stats=True,
                         debug=False):

    """
    create a figure comparing 1D burial and thermal model results
    with vitrinite reflectance, apatite fission track and present-day
    temperature data

    :param model_run_data:
    :param show_provenance_hist:
    :param show_strat_column:
    :param show_thermochron_data:
    :param contour_variable:
    :param add_legend:
    :param strat_fontsize:
    :param figsize:
    :param legend_space:
    :param height_ratio:
    :param cb_buffer_vert:
    :param cb_buffer_hor:
    :param show_prov_ages_simple:
    :param ncols_legend:
    :param add_panel_titles: add a label to each figure panel
    :param panel_title_level: if 1 label panels a, b, c. If 2, label panels as a1, a2, a3 etc...
    :param panel_label_fs: fontsize of panel label
    :param contour_variable: 'temperature' (default), 'salinity' or
        'pressure'. selects which variable is shown in the main burial
        history panel; 'pressure' shows modeled excess (above
        hydrostatic) pore pressure, only available for model runs with
        simulate_fluid_flow enabled
    :param show_porosity_panel: add a panel with modeled (and, if
        available, observed) present-day porosity vs depth. defaults
        (None) to True if this model run has excess pressure data
        (simulate_fluid_flow was used), False otherwise
    :param show_pressure_panel: add a panel with modeled (and, if
        available, observed) present-day pressure vs depth. only shown
        if the model run has excess pressure data (simulate_fluid_flow).
        defaults (None) to True if this model run has excess pressure
        data, False otherwise
    :param show_gof_stats: choose which model fit statistics (GOF, RMSE
        and/or error) are added to the figure panels. ``True`` (default)
        shows the fit statistics for every data panel that is present,
        ``False`` hides all of them. Pass a list or tuple of panel names
        instead to select individual panels, choosing from 'temperature',
        'VR', 'AFT', 'AHe' and 'pressure'.
    :param debug:
    :return:
    """

    # accept either the current named dict/xarray output format or the
    # older plain positional list format (older .pck files); normalize
    # both to the same named format before unpacking
    normalized = pybasin_io.normalize_model_run_data(model_run_data)

    grid = normalized['grid']
    time_array_bp = grid['time'].values
    surface_temp_array = grid['surface_temperature'].values
    basal_hf_array = grid['basal_heat_flow'].values
    z_nodes = grid['z'].values
    active_nodes = grid['active'].values
    T_nodes = grid['T'].values
    node_strat = grid['node_strat'].values.tolist()
    node_age = grid['node_age'].values

    # porosity is always calculated, but older saved model runs may not
    # have it; excess pore pressure is only calculated for model runs
    # with simulate_fluid_flow enabled
    porosity_nodes = grid['porosity'].values if 'porosity' in grid else None
    rho_nodes = grid['rho'].values if 'rho' in grid else None
    P_ex_nodes = grid['P_ex'].values if 'P_ex' in grid else None

    pressure_data = normalized.get('pressure_data')
    if pressure_data is not None:
        # older saved model runs only carried the raw observed-data
        # DataFrame for overlay, with no goodness-of-fit statistics
        pressure_gof = pressure_data.get('gof', np.nan)
        pressure_rmse = pressure_data.get('rmse', np.nan)
    porosity_data = normalized.get('porosity_data')

    # by default, show the porosity and pressure panels once
    # compaction-driven pore pressure modeling (simulate_fluid_flow)
    # was used for this model run, since the excess pressure and
    # present-day pressure comparison are then directly relevant
    if show_porosity_panel is None:
        show_porosity_panel = P_ex_nodes is not None
    if show_pressure_panel is None:
        show_pressure_panel = P_ex_nodes is not None

    if show_pressure_panel and P_ex_nodes is None:
        logger.info('show_pressure_panel is True, but this model run has no '
                   'excess pressure data (simulate_fluid_flow was not used); '
                   'not showing the pressure panel')
        show_pressure_panel = False

    if contour_variable == 'pressure' and P_ex_nodes is None:
        logger.info('contour_variable is "pressure", but this model run has '
                   'no excess pressure data (simulate_fluid_flow was not '
                   'used); showing temperature instead')
        contour_variable = 'temperature'

    # normalize show_gof_stats to a set of panel names for which fit
    # statistics should be shown: True/False select all/none, a list or
    # tuple selects individual panels ('temperature', 'VR', 'AFT',
    # 'AHe', 'pressure')
    if show_gof_stats is True:
        show_gof_stats = {'temperature', 'VR', 'AFT', 'AHe', 'pressure'}
    elif show_gof_stats is False:
        show_gof_stats = set()
    else:
        show_gof_stats = set(show_gof_stats)

    T_data = normalized['T_data']
    if T_data is not None:
        T_depth = T_data['depth']
        T_obs = T_data['observed']
        T_obs_sigma = T_data['observed_sigma']
        T_data_type = T_data['data_type']
        T_gof = T_data['gof']
        T_rmse = T_data['rmse']

    C_data = normalized['C_data']
    if C_data is not None:
        C_nodes = C_data['C_nodes']
        surface_salinity_array = C_data['surface_salinity']
        salinity_lwr_bnd = C_data['salinity_lower_bnd']
        salinity_depth = C_data['salinity_depth']
        salinity_data = C_data['salinity_observed']
        salinity_data_unc = C_data['salinity_observed_unc']
        salinity_RMSE = C_data['salinity_rmse']
        q_solute_bottom = C_data['q_solute_bottom']
        q_solute_top = C_data['q_solute_top']

    VR_model_data = normalized['VR_model_data']
    if VR_model_data is not None:
        vr_nodes = VR_model_data['vr_nodes']
        vr_depth = VR_model_data['depth']
        vr_obs = VR_model_data['observed']
        vr_min = VR_model_data['observed_min']
        vr_max = VR_model_data['observed_max']
        vr_obs_sigma = VR_model_data['observed_sigma']
        vr_GOF = VR_model_data['gof']
        vr_rmse = VR_model_data['rmse']
        vr_data_well = VR_model_data['data']

    AFT_data = normalized['AFT_data']
    if AFT_data is not None:
        simulated_AFT_data = AFT_data['simulated']
        aft_sample_names = AFT_data['sample_names']
        aft_age_depth = AFT_data['age_depth']
        aft_age = AFT_data['age']
        aft_age_stderr_min = AFT_data['age_stderr_min']
        aft_age_stderr_plus = AFT_data['age_stderr_plus']
        aft_length_mean = AFT_data['length_mean']
        aft_length_std = AFT_data['length_std']
        aft_age_samples = AFT_data['age_samples']
        single_grain_aft_ages = AFT_data['single_grain_ages']
        single_grain_aft_ages_se_min = AFT_data['single_grain_ages_se_min']
        single_grain_aft_ages_se_plus = AFT_data['single_grain_ages_se_plus']
        aft_age_bins = AFT_data['age_bins']
        aft_age_pdfs = AFT_data['age_pdfs']
        aft_age_GOF = AFT_data['gof']
        aft_age_error = AFT_data['age_error']
        aft_sample_times = AFT_data['sample_times']
        aft_sample_temps = AFT_data['sample_temps']
        # note: time_array_bp is intentionally overwritten here, matching
        # the historical positional format where AFT_data carried its own copy
        time_array_bp = AFT_data['time_array_bp']
        z_aft_samples = AFT_data['z_samples']
        T_samples = AFT_data['T_samples']
        aft_data_samples = AFT_data['data']

    AHe_data = normalized['AHe_data']
    if AHe_data is not None:
        ahe_sample_depths = AHe_data['sample_depths']
        ahe_ages_all_samples = AHe_data['ages_all_samples']
        ahe_ages_all_samples_SE = AHe_data['ages_all_samples_se']
        ahe_age_bin = AHe_data['age_bin']
        ahe_age_pdfs = AHe_data['age_pdfs']
        modeled_ahe_age_samples = AHe_data['modeled_age_samples']
        modeled_ahe_age_samples_min = AHe_data['modeled_age_samples_min']
        modeled_ahe_age_samples_max = AHe_data['modeled_age_samples_max']
        ahe_age_gof = AHe_data['gof']
        ahe_age_error = AHe_data['age_error']
        simulated_AHe_data = AHe_data['simulated']
        ahe_data_samples = AHe_data['data']

    nt_total, n_nodes = T_nodes.shape

    if AFT_data is not None and simulated_AFT_data is not None:
        aft_age_nodes = simulated_AFT_data['age_nodes']
        aft_age_nodes_min = simulated_AFT_data['age_nodes_min']
        aft_age_nodes_max = simulated_AFT_data['age_nodes_max']
        aft_ln_mean_nodes = simulated_AFT_data['ln_mean_nodes']
        aft_ln_std_nodes = simulated_AFT_data['ln_std_nodes']
        aft_node_times_burial = simulated_AFT_data['node_times_burial']
        aft_node_zs = simulated_AFT_data['node_zs']
        aft_node_times = simulated_AFT_data['node_times']
        aft_node_temps = simulated_AFT_data['node_temps']

        _, n_prov_scenarios, n_kinetic_scenarios = aft_age_nodes.shape

        prov_ages = [aft_node_times[0][0].max(),
                     aft_node_times[0][1].max()]


    if AHe_data is not None and simulated_AHe_data is not None:
        ahe_age_nodes = simulated_AHe_data['age_nodes']
        ahe_age_nodes_min = simulated_AHe_data['age_nodes_min']
        ahe_age_nodes_max = simulated_AHe_data['age_nodes_max']
        ahe_node_times_burial = simulated_AHe_data['node_times_burial']
        ahe_node_zs = simulated_AHe_data['node_zs']

        prov_ages = [ahe_node_times_burial[0][0].max(),
                     ahe_node_times_burial[0][-1].max()]

    PY3 = sys.version_info.major == 3

    if PY3:
        degree_symbol = chr(176)
    else:
        degree_symbol = unichr(176)

    xsize = figsize / 25.4
    golden_ratio = (1.0 + np.sqrt(5))/2.0
    ysize = xsize / golden_ratio

    if add_legend is True:
        ysize = ysize * 1.12

    font = {'family': 'sans-serif',
            'size': 9}

    pl.rc('font', **font)

    width_ratios = [8]

    ncols = 1

    if add_legend is True:
        bottom += legend_space

    leg_items = []
    leg_labels = []
    model_label = []
    model_range_label = []
    data_label = []
    data_ext_label = []
    leg_data = None
    leg_data_ext = []
    leg_model_range = None

    max_depth = z_nodes[active_nodes].max() * 1.1

    # skip VR, AFT and AHe panels if no data
    #if VR_model_data is not None and len(vr_obs) == 0:
    #    VR_model_data = None
    #if AFT_data is not None and len(aft_age) == 0:
    #    AFT_data = None
    #if AHe_data is not None and len(ahe_ages_all_samples) == 0:
    #    AHe_data = None

    if show_thermochron_data is False:
        logger.info('not showing thermochron data:')
        AFT_data = None
        AHe_data = None

    aft_sg_sorted = []
    if show_age_scatter and AFT_data is not None:
        aft_sg_sorted = sorted(
            [j for j in range(len(aft_age_depth))
             if single_grain_aft_ages[j] is not None
             and len(single_grain_aft_ages[j]) > 0],
            key=lambda j: aft_age_depth[j],
        )

    ahe_sorted_indices = []
    if show_age_scatter and AHe_data is not None and modeled_ahe_age_samples is not None:
        ahe_sorted_indices = sorted(
            range(len(ahe_sample_depths)),
            key=lambda i: ahe_sample_depths[i],
        )

    n_aft_sc = len(aft_sg_sorted)
    n_ahe_sc = len(ahe_sorted_indices)

    if n_aft_sc > 0 or n_ahe_sc > 0:
        # a shared, fine-grained row grid lets the stacked scatter panels
        # each get an equal, whole number of rows without nesting a
        # subfigure or sub-gridspec for them: nesting shrinks axes that
        # use set_box_aspect to a fraction of their allotted space under
        # a constrained-layout figure
        if n_aft_sc > 0 and n_ahe_sc > 0:
            n_scatter_units = math.lcm(n_aft_sc, n_ahe_sc)
        else:
            n_scatter_units = n_aft_sc or n_ahe_sc

        # extra resolution so the gap reserved between stacked scatter
        # panels (for their title and tick labels) can be tuned finely.
        # constrained_layout solve time grows fast with the row count, so
        # keep this modest
        row_resolution = 2
        n_scatter_units *= row_resolution

        hr_units = max(int(round(height_ratio)), 1)
        row_units = 1 + hr_units + 1
        nrows = row_units * n_scatter_units
        row_top = slice(0, n_scatter_units)
        row_mid = slice(n_scatter_units, n_scatter_units * (1 + hr_units))
        row_bot = slice(n_scatter_units * (1 + hr_units), nrows)
        row_height_ratios = [1] * nrows
    else:
        nrows = 3
        row_top, row_mid, row_bot = 0, 1, 2
        row_height_ratios = [1, height_ratio, 1]

    temp_panel_ind = None
    if show_strat_column is True:
        strat_panel_ind = ncols
        ncols += 1
        width_ratios.append(2)
        if show_temp_panel:
            temp_panel_ind = ncols
            ncols += 1
            width_ratios.append(3)
    elif show_temp_panel:
        temp_panel_ind = ncols
        ncols += 1
        width_ratios.append(3)

    porosity_panel_ind = None
    if show_porosity_panel:
        porosity_panel_ind = ncols
        ncols += 1
        width_ratios.append(3)

    pressure_panel_ind = None
    if show_pressure_panel:
        pressure_panel_ind = ncols
        ncols += 1
        width_ratios.append(3)

    if C_data is not None:
        C_panel_ind = ncols
        ncols += 1
        width_ratios.append(3)

    if VR_model_data is not None:
        logger.info('adding panel for VR data')
        vr_panel_ind = ncols
        ncols += 1
        width_ratios.append(3)

    if AFT_data is not None:
        logger.info('adding panel for AFT data')
        aft_panel_ind = ncols
        ncols += 1
        width_ratios.append(4)

    if AHe_data is not None:
        logger.info('adding panel for AHe data')
        ahe_panel_ind = ncols
        ncols += 1
        width_ratios.append(4)

    # width of the panels that the default figsize was tuned for, i.e.
    # everything up to and excluding the age scatter columns
    core_width_units = sum(width_ratios)

    aft_scatter_col = None
    ahe_scatter_col = None
    if show_age_scatter:
        has_scatter = len(aft_sg_sorted) > 0 or len(ahe_sorted_indices) > 0
        if has_scatter:
            # spacer column so y-axis labels of scatter panels do not overlap the age panels
            ncols += 1
            width_ratios.append(2)
        if len(aft_sg_sorted) > 0:
            aft_scatter_col = ncols
            ncols += 1
            width_ratios.append(4)
        if len(ahe_sorted_indices) > 0:
            ahe_scatter_col = ncols
            ncols += 1
            width_ratios.append(4)
        if has_scatter:
            # trailing spacer so the x-axis label of the rightmost scatter
            # panel is not clipped at the edge of the figure
            ncols += 1
            width_ratios.append(1)

    # widen the figure to fit the extra scatter columns
    xsize = xsize * sum(width_ratios) / core_width_units

    # each stacked scatter panel needs a minimum height to remain
    # readable, so grow the figure height for wells with several samples
    n_scatter_rows = max(
        len(aft_sg_sorted) if aft_scatter_col is not None else 0,
        len(ahe_sorted_indices) if ahe_scatter_col is not None else 0,
    )
    if n_scatter_rows > 0:
        min_scatter_height = min(n_scatter_rows * 1.6, 14.0)
        ysize = max(ysize, min_scatter_height)

    fig = pl.figure(figsize=(xsize, ysize), layout="constrained")

    gs = gridspec.GridSpec(nrows, ncols,
                           wspace=0.06, hspace=0.08,
                           bottom=bottom, top=top,
                           left=left, right=right,
                           width_ratios=width_ratios,
                           height_ratios=row_height_ratios,
                           figure=fig)

    axb = fig.add_subplot(gs[row_mid, 0])
    axst = fig.add_subplot(gs[row_top, 0])
    axhf = fig.add_subplot(gs[row_bot, 0])

    all_panels = [axst, axb, axhf]

    if show_strat_column is True:
        ax_strat = fig.add_subplot(gs[row_mid, strat_panel_ind])
        all_panels.append(ax_strat)

    if show_temp_panel:
        ax_temp = fig.add_subplot(gs[row_mid, temp_panel_ind])
        all_panels.append(ax_temp)

    if show_porosity_panel:
        ax_porosity = fig.add_subplot(gs[row_mid, porosity_panel_ind])
        all_panels.append(ax_porosity)

    if show_pressure_panel:
        ax_pressure = fig.add_subplot(gs[row_mid, pressure_panel_ind])
        all_panels.append(ax_pressure)

    if C_data is not None:
        ax_c = fig.add_subplot(gs[row_mid, C_panel_ind])
        all_panels.append(ax_c)
    if VR_model_data is not None:
        ax_vr = fig.add_subplot(gs[row_mid, vr_panel_ind])
        all_panels.append(ax_vr)
    if AFT_data is not None:
        ax_afta = fig.add_subplot(gs[row_mid, aft_panel_ind])
        all_panels.append(ax_afta)
    if AHe_data is not None:
        ax_ahe = fig.add_subplot(gs[row_mid, ahe_panel_ind])
        all_panels.append(ax_ahe)

    def scatter_row_slices(n_samples, rows_per_panel):
        # reserve part of each panel's rows as a gap for the title and
        # tick labels of the next panel down
        gap = max(int(round(rows_per_panel * 0.3)), 1) if rows_per_panel > 1 else 0
        return [slice(i * rows_per_panel, (i + 1) * rows_per_panel - gap)
                for i in range(n_samples)]

    # scatter panels are placed directly on the shared, fine-grained gs
    # rather than in a nested subfigure or sub-gridspec, since nesting
    # shrinks axes that use set_box_aspect under constrained layout
    ax_aft_scatter_list = []
    ax_ahe_scatter_list = []
    if show_age_scatter:
        if aft_scatter_col is not None:
            rows_per_aft = nrows // n_aft_sc
            aft_row_slices = scatter_row_slices(n_aft_sc, rows_per_aft)
            ax_aft_scatter_list = [fig.add_subplot(gs[s, aft_scatter_col]) for s in aft_row_slices]
        if ahe_scatter_col is not None:
            rows_per_ahe = nrows // n_ahe_sc
            ahe_row_slices = scatter_row_slices(n_ahe_sc, rows_per_ahe)
            ax_ahe_scatter_list = [fig.add_subplot(gs[s, ahe_scatter_col]) for s in ahe_row_slices]

    depth_panels = [all_panels[1]] + all_panels[3:]
    time_panels = all_panels[:3]

    line_props = {"color": "black", "lw": 1.0}

    scatter_props = {"marker": "o",
                     "s": 25,
                     "color": 'gray',
                     "edgecolor": 'black',
                     "zorder": 10}

    erb_props = {"marker": "o",
                 "ms": 4,
                 "linestyle": "None",
                 "color": 'black',
                 "mec": 'black',
                 "mfc": 'gray',
                 "lw": 0.75,
                 "zorder": 10}

    textprops = {"fontsize": 'small',
                 'ha': 'center',
                 'va': 'bottom',
                 'weight': 'normal',
                 'bbox': dict(facecolor="white",
                              ec='white',
                              alpha=0.7)}

    provenance_color = 'darkgray'
    cmap = matplotlib.cm.get_cmap('coolwarm')

    if contour_variable == 'salinity':
        cnt_var = C_nodes
        cnt_step = 0.005
        cb_label = 'salinity (kg/kg)'
    elif contour_variable == 'pressure':
        cnt_var = P_ex_nodes / 1.0e6
        cnt_step = max(cnt_var[active_nodes].max() / 15.0, 0.1)
        cb_label = 'excess pressure (MPa)'
    else:
        cnt_var = T_nodes
        if T_nodes.max() < 50:
            cnt_step = 2.5
        elif T_nodes.max() < 100:
            cnt_step = 5.0
        else:
            cnt_step = 10.0
        cb_label = 'T (%s C)' % degree_symbol

    # plot surface temperature. no equivalent surface line for
    # pressure: excess pressure is 0 at the surface by construction
    if contour_variable == 'salinity':
        axst.plot(time_array_bp / 1e6, surface_salinity_array,
                  **line_props)
    elif contour_variable != 'pressure':
        axst.plot(time_array_bp / 1e6, surface_temp_array,
                  **line_props)

    ts = 1.0e5

    if max_depth < 1000:
        ys = 1.0
    else:
        ys = 10.0

    yi = np.arange(z_nodes[active_nodes].min(), max_depth + ys, ys)

    cnt_var_mask = cnt_var.copy()
    cnt_var_mask[active_nodes == False] = np.nan

    time_2d = np.zeros([nt_total, n_nodes])

    for nn in range(n_nodes):
        time_2d[:, nn] = time_array_bp / 1.0e6

    #
    mean_timestep = np.mean(-np.diff(time_array_bp))
    time_int_grid = int(np.round(ts / mean_timestep))

    ntsx = len(time_2d[::time_int_grid])

    xi = np.linspace(np.min(time_array_bp), np.max(time_array_bp), ntsx) \
        / 1.0e6

    x = time_2d[::time_int_grid].ravel()
    y = z_nodes[::time_int_grid].ravel()
    z = cnt_var_mask[::time_int_grid].ravel()
    act = active_nodes[::time_int_grid].ravel()
    ind_act = act == True

    logger.info('gridding T or salinity data vs time')
    gridding_ok = True
    # serial 1D interpolation, failproof method, 2D interpolation fails or
    # inaccurate with strongly different x,y scales
    zi_data = np.zeros((len(yi), len(xi)))
    zi = np.ma.masked_array(zi_data, mask=np.isnan(zi_data))
    nts = len(time_2d[::time_int_grid])
    for tsi in range(nts):
        y_1d = z_nodes[tsi * time_int_grid]
        z_1d = cnt_var_mask[tsi * time_int_grid]
        ind_nan = np.isnan(z_1d) == False
        z_interpolated = np.interp(yi, y_1d[ind_nan], z_1d[ind_nan])
        zi[:, -tsi] = z_interpolated

    if gridding_ok is True:
        # find max depth at each timestep
        z_nodes_corr = z_nodes.copy()
        z_nodes_corr[np.isnan(z_nodes_corr)] = -99999
        max_depth_time = np.max(z_nodes_corr, axis=1)
        max_depth_time2 = np.interp(xi, (time_array_bp/1.0e6)[::-1], max_depth_time[::-1])

        # filter interpolated values that are deeper than deepest fm.
        for nti in range(len(xi)):
            zi.mask[yi > max_depth_time2[nti], nti] = True

        #tc = axb.pcolormesh(xg, yg, zi2, cmap='jet')
        c_int = np.arange(0.0, cnt_var[active_nodes].max()+cnt_step, cnt_step)
        tc = axb.contourf(xi, yi, zi, c_int, cmap=cmap, zorder=1.0)

    else:
        plot_int = 1
        tc = axb.scatter(x[ind_act][::plot_int],
                         y[ind_act][::plot_int],
                         c=z[ind_act][::plot_int],
                         edgecolor="black", lw=0.1,
                         s=10,
                         cmap=cmap)
        
    logger.debug(f"node strat: {node_strat}")

    major_strat = [n.split('_s_')[0].split('_a_')[0] for n in node_strat]

    #major_strat = [msi for msi in major_strat]

    #major_strat = node_strat.copy()    
    strat_transition = [m != n for m, n in zip(major_strat[:-1],
                                               major_strat[1:])]
    strat_transition.append(True)

    strat_transition = np.array(strat_transition)

    # check and reduce number of strat units shown
    n_strat_units_shown = np.sum(strat_transition)

    if max_strat_units is not None and n_strat_units_shown > max_strat_units:
        logger.info('reducing number of strat units shown from %i to %i' % (n_strat_units_shown, max_strat_units))
        sint = int(np.ceil(n_strat_units_shown / max_strat_units))

        ind = np.where(strat_transition == True)[0]
        strat_transition[:] = False
        strat_transition[ind[::sint]] = True
        strat_transition[ind[-1]] = True

    units_shown = [major_strat[i] for i, s in enumerate(strat_transition) if s == True]
    logger.debug('strat units shown in fig: %s' % units_shown)

    # plot provenance and burial histories
    #if (AFT_data is not None or AHe_data is not None) \
    #        and show_provenance_hist is True:
    if ((AFT_data is not None and simulated_AFT_data is not None) or (AHe_data is not None and simulated_AHe_data is not None)) \
            and show_provenance_hist is True:

        # get burial histories for all the thermochron nodes

        if AFT_data is not None:
            burial = aft_node_times_burial
            depths = aft_node_zs
        else:
            burial = ahe_node_times_burial
            depths = ahe_node_zs

        # show provenance and burial histories
        strat_count = 0
        for xb, yb, strat_trans in zip(burial, depths,
                                       strat_transition):
            if strat_trans == True:
                c = provenance_color
                cf = 'beige'

                # find min and max provenance depth
                min_prov_ind = 0
                max_prov_ind = 0
                max_prov_age = 0
                min_prov_age = 99999
                for i, xbi, ybi in zip(itertools.count(), xb, yb):
                    ind_surface = np.where(ybi >= 0)[0][0]
                    prov_age_mid = np.mean(xbi[ind_surface])
                    if prov_age_mid > max_prov_age:
                        max_prov_age = prov_age_mid
                        max_prov_ind = i
                    if prov_age_mid < min_prov_age:
                        min_prov_age = prov_age_mid
                        min_prov_ind = i

                #xf = np.concatenate((xb[min_prov_ind], xb[max_prov_ind][::-1]))
                #yf = np.concatenate((yb[min_prov_ind], yb[max_prov_ind][::-1]))

                #axb.fill(xf, yf, color=cf, zorder=0)
                #leg_prov_fill = mpatches.Patch(color=cf)

                combs = list(itertools.combinations(list(range(len(xb))), 2))
                for comb in combs:

                    i1, i2 = comb
                    xf = np.concatenate((xb[i1], xb[i2][::-1]))
                    yf = np.concatenate((yb[i1], yb[i2][::-1]))

                    axb.fill(xf, yf, color=cf, zorder=0)
                    leg_prov_fill = mpatches.Patch(color=cf)

                for xbi, ybi in zip(xb, yb):
                    leg_prov, = axb.plot(xbi, ybi, color=c, lw=0.5)

                strat_count += 1

        leg_items += [leg_prov, leg_prov_fill]
        leg_labels += ['provenance and burial history',
                       'range of provenance histories']

    else:
        ind = np.array(strat_transition) == True
        n_strat_trans = ind.sum()
        for i in range(n_strat_trans):
            leg_strat_unit, = axb.plot(time_array_bp / 1e6,
                                       z_nodes[:, ind][:, i],
                                       color='black',
                                       lw=0.5, zorder=100)

        leg_items += [leg_strat_unit]
        leg_labels += ['stratigraphic unit']

    if (AFT_data is not None or AHe_data is not None or show_thermochron_data is False) \
            and show_prov_ages_simple is True:

        logger.info('showing errorbar for AFT start times:')

        x = np.array(prov_ages).mean()
        xerr = np.abs(x - prov_ages[0])
        #leg_prov_simple = axb.scatter(prov_ages, [0, 0],
        #                              marker='*',
        #                              facecolor='gray', edgecolor='black',
        #                              zorder=301)
        leg_prov_simple = axb.errorbar([x], [0], xerr=[xerr],
                                       marker='None', color='gray',
                                       lw=1.0,
                                       zorder=301)
        leg_items += [leg_prov_simple]
        leg_labels += ['provenance ages']

    # plot basal heat flow. no equivalent lower boundary line for
    # pressure (a no-flow boundary), so nothing is plotted here
    if contour_variable == 'salinity':
        axhf.axhline(y=salinity_lwr_bnd, **line_props)
    elif contour_variable != 'pressure':
        axhf.plot(time_array_bp / 1e6, basal_hf_array * 1000.0,
                  **line_props)
        axhf.set_ylim(basal_hf_array.min() * 1000.0 * 0.95,
                      basal_hf_array.max() * 1000.0 * 1.05)

    if show_temp_panel:
        leg_model, = ax_temp.plot(T_nodes[-1, active_nodes[-1]],
                                  z_nodes[-1, active_nodes[-1]],
                                  **line_props)
        model_label.append('temperature')

    if show_porosity_panel:
        ax_porosity.plot(porosity_nodes[-1, active_nodes[-1]],
                         z_nodes[-1, active_nodes[-1]],
                         **line_props)
        model_label.append('porosity')

        if porosity_data is not None:
            obs = porosity_data['data']
            ax_porosity.scatter(obs['porosity'], obs['depth'], **scatter_props)

    if show_pressure_panel:
        # hydrostatic pore water pressure (1025 kg/m3, 9.81 m/s2),
        # matching the assumption used throughout the compaction and
        # fluid flow model
        z_active = z_nodes[-1, active_nodes[-1]]
        hydrostatic_pressure = 1025.0 * 9.81 * z_active / 1.0e6
        total_pressure = hydrostatic_pressure \
            + P_ex_nodes[-1, active_nodes[-1]] / 1.0e6

        ax_pressure.plot(hydrostatic_pressure, z_active,
                         color='gray', lw=1.0, ls='--',
                         label='hydrostatic')

        # lithostatic (total overburden) pressure, from the bulk
        # density already calculated for the heat flow model. only
        # available for model runs that saved bulk density (rho);
        # older saved model runs may not have it
        lithostatic_pressure = None
        if rho_nodes is not None:
            rho_active = rho_nodes[-1, active_nodes[-1]]
            lithostatic_pressure = np.zeros_like(z_active)
            lithostatic_pressure[1:] = np.cumsum(
                0.5 * (rho_active[:-1] + rho_active[1:]) * 9.81
                * np.diff(z_active)) / 1.0e6
            ax_pressure.plot(lithostatic_pressure, z_active,
                             color='gray', lw=1.0, ls=':',
                             label='lithostatic')

        ax_pressure.plot(total_pressure, z_active, **line_props)
        model_label.append('pressure')

        if pressure_data is not None:
            obs = pressure_data['data']
            obs_depth = (obs['depth_top'] + obs['depth_bottom']) / 2.0
            ax_pressure.scatter(obs['FSIP_MPa'], obs_depth, **scatter_props)

    if show_strat_column is True:

        ind = np.array(strat_transition) == True
        n_strat_trans = ind.sum()

        z_trans = z_nodes[:, ind][-1, :]
        z_trans = np.insert(z_trans, [0], np.array([0]))
        strat_trans = np.array(major_strat)[ind]
        for ax in depth_panels[1:]:
            for i in range(n_strat_trans):
                leg_strat_unit = ax.axhline(y=z_nodes[:, ind][-1, i],
                                            color='gray',
                                            lw=0.5, zorder=1)

        # add labels for stratigraphic units
        z_mid_trans = (z_trans[1:] + z_trans[:-1]) / 2.0
        for z_pos, strat_name in zip(z_mid_trans, strat_trans):
            if strat_name[0] != "+":
                ax_strat.text(0.03, z_pos, strat_name, fontsize=strat_fontsize)

    if show_temp_panel and T_data is not None and len(T_data) > 0:
        ind = T_data_type == 'BHT'
        nind = T_data_type != 'BHT'

        ind = ind.values
        nind = nind.values

        if 'BHT' in T_data_type.values:
            xerr = np.array([np.zeros_like(T_obs_sigma)[ind], T_obs_sigma[ind] * 2])
            leg_data = ax_temp.errorbar(T_obs[ind], T_depth[ind], xerr=xerr, **erb_props)

        leg_data = ax_temp.errorbar(T_obs[nind], T_depth[nind], xerr=T_obs_sigma[nind] * 2, **erb_props)
        data_label.append('temperature')

    # plot modeled salinity
    if C_data is not None and C_nodes is not None:
        leg_model, = ax_c.plot(C_nodes[-1, active_nodes[-1]],
                               z_nodes[-1, active_nodes[-1]],
                               **line_props)
        model_label.append('salinity')

    if C_data is not None and len(salinity_data) > 0:
        leg_data = ax_c.scatter(salinity_data, salinity_depth,
                                **scatter_props)
        data_label.append('salinity')

    # plot vitrinite
    if VR_model_data is not None and vr_nodes is not None:
        leg_model, = ax_vr.plot(vr_nodes[-1, active_nodes[-1]],
                                z_nodes[-1, active_nodes[-1]],
                                **line_props)

    if VR_model_data is not None and len(vr_obs) > 0:
        #if VR_model_data is not None and len(vr_data) > 0:
        xerr = np.ones((2, len(vr_depth))) * vr_obs_sigma

        ind = np.isnan(vr_min) == False
        xerr[0][ind] = vr_obs[ind] - vr_min[ind]
        xerr[1][ind] = vr_max[ind] - vr_obs[ind]

        leg_data = ax_vr.errorbar(vr_obs, vr_depth,
                                  xerr=xerr,
                                  **erb_props)

        model_label.append('VR')
        data_label.append('VR')

    # plot modeled aft ages
    if AFT_data is not None and simulated_AFT_data is not None:
        ax_afta.fill_betweenx(z_nodes[-1, active_nodes[-1]],
                              aft_age_nodes_min[active_nodes[-1]],
                              aft_age_nodes_max[active_nodes[-1]],
                              color='lightgrey')
        leg_model_range = mpatches.Patch(color='lightgrey')

        leg_strat, = ax_afta.plot(node_age[active_nodes[-1]],
                                  z_nodes[-1, active_nodes[-1]],
                                  color='green', lw=1.5, ls='--', zorder=101)
        leg_items.append(leg_strat)
        leg_labels.append('age of deposition')

        model_range_label.append('AFT ages')

    if AFT_data is not None:

        # violin plots of single grain age pdf
        if show_violin_plot is True:
            violin_width = max_depth / 20.0
            pdf_threshold = 1e-5
            for sample_no in range(len(aft_age)):
                pdf_plot = aft_age_pdfs[sample_no]

                if np.any(np.isnan(pdf_plot)) == False and \
                                single_grain_aft_ages[sample_no] is not None:

                    ind = pdf_plot > pdf_threshold

                    #pdf_plot[pdf_plot < pdf_threshold] = 0.0
                    vd = dict(coords=aft_age_bins[sample_no][ind],
                              vals=aft_age_pdfs[sample_no][ind],
                              mean=1.0, min=1.0, max=1.0, median=1.0)
                    vp = ax_afta.violin([vd],
                                        positions=[aft_age_depth[sample_no]],
                                        vert=False,
                                        widths=violin_width,
                                        showextrema=False)
                    for pc in vp['bodies']:
                        pc.set_edgecolor('black')
                        pc.set_facecolor('lightblue')
                        pc.set_alpha(0.75)
                        pc.set_linewidth(0.5)

            leg_violin = mpatches.Patch(facecolor='lightblue',
                                        edgecolor='black', lw=0.5)
            leg_data_ext.append(leg_violin)
            data_ext_label.append('age distribution')

        # show single grain AFT ages, without errorbar
        for sample_no in range(len(single_grain_aft_ages)):
            x = single_grain_aft_ages[sample_no]

            if x is not None:
                y = np.ones_like(x) * aft_age_depth[sample_no]
                leg_sg = ax_afta.scatter(x, y, color='black', s=5, marker='o')
                if len(data_ext_label) == 0 or 'single grain AFT ages' not in data_ext_label:
                    leg_data_ext.append(leg_sg)
                    data_ext_label.append('single grain AFT ages')

        #ind_ca = np.array([a is None for a in single_grain_aft_ages])
        ind_ca = np.array([True] * len(single_grain_aft_ages))

        if True in ind_ca:
            # show central ages
            leg_data = ax_afta.errorbar(aft_age[ind_ca], aft_age_depth[ind_ca],
                                        xerr=[aft_age_stderr_min[ind_ca] * 1.96,
                                              aft_age_stderr_plus[ind_ca] * 1.96],
                                        **erb_props)
            #if len(leg_labels) == 0 or 'AFT age' not in leg_labels[-1]:
            data_label.append('AFT age')

    if AFT_data is not None and simulated_AFT_data is not None:
        for n_prov in range(n_prov_scenarios):
            for n_kin in range(n_kinetic_scenarios):
                leg_model, = ax_afta.plot(
                    aft_age_nodes[active_nodes[-1], n_prov, n_kin],
                    z_nodes[-1, active_nodes[-1]],
                    **line_props)
        model_label.append('AFT ages')

    # plot track lengths
    #for n_prov in xrange(n_prov_scenarios):
    #    for n_kin in xrange(n_kinetic_scenarios):
    #        ax_aftln.fill_betweenx(z_nodes[-1, active_nodes[-1]],
    #                               (aft_ln_mean_nodes[active_nodes[-1], n_prov, n_kin]
    #                                - aft_ln_std_nodes[active_nodes[-1], n_prov, n_kin]),
    #                               (aft_ln_mean_nodes[active_nodes[-1], n_prov, n_kin]
    #                                + aft_ln_std_nodes[active_nodes[-1], n_prov, n_kin]),
    #                               color='lightgrey', zorder=0)

    #        ax_aftln.plot(aft_ln_mean_nodes[active_nodes[-1], n_prov, n_kin],
    #                      z_nodes[-1, active_nodes[-1]], zorder=10,
    #                      **line_props)

    #ax_aftln.errorbar(aft_length_mean, aft_age_depth,
    #                  xerr=aft_length_std, **erb_props)

    if AHe_data is not None and simulated_AHe_data is not None:

        _, n_grain_radius, n_prov_scenarios = np.array(ahe_age_nodes).shape

        ahe_age_min_grains = np.array(ahe_age_nodes_min)[active_nodes[-1]]
        ahe_age_max_grains = np.array(ahe_age_nodes_max)[active_nodes[-1]]

        ahe_age_min = np.min(ahe_age_min_grains, axis=1)
        ahe_age_max = np.max(ahe_age_max_grains, axis=1)

        ax_ahe.fill_betweenx(z_nodes[-1, active_nodes[-1]],
                             ahe_age_min, ahe_age_max,
                             color='lightgrey')
        leg_model_range = mpatches.Patch(color='lightgrey')
        leg_strat, = ax_ahe.plot(node_age[active_nodes[-1]],
                                 z_nodes[-1, active_nodes[-1]],
                                 color='green', lw=1.5, ls='--', zorder=101)
        leg_data_ext.append(leg_strat)
        data_ext_label.append('age of deposition')

        ahe_age_nodes_array = np.array(ahe_age_nodes)
        for n_prov in range(n_prov_scenarios):
            for n_rad in range(n_grain_radius):
                leg_model, = ax_ahe.plot(ahe_age_nodes_array[active_nodes[-1], n_rad, n_prov],
                                         z_nodes[-1, active_nodes[-1]],
                                         **line_props)

        model_label.append('He ages')
        model_range_label.append('He ages')

    # show modelled age envelope for samples
    if AHe_data is not None and modeled_ahe_age_samples is not None:

        modelled_ahe_age_min_all_samples = np.min(np.array(modeled_ahe_age_samples_min), axis=1)
        modelled_ahe_age_max_all_samples = np.max(np.array(modeled_ahe_age_samples_max), axis=1)
        depths_array= np.array(ahe_sample_depths)

        leg_ahe_sample_min = ax_ahe.scatter(modelled_ahe_age_min_all_samples, depths_array, marker="|", zorder=100, color="tab:blue")
        leg_ahe_sample_max = ax_ahe.scatter(modelled_ahe_age_max_all_samples, depths_array, marker="|", zorder=100, color="tab:blue")
    else:
        leg_ahe_sample_min = None
        leg_ahe_sample_max = None
        
    # show AHe data
    if AHe_data is not None:
        for ahe_ages_sample, ahe_sample_depth, ahe_ages_sample_SE in \
                zip(ahe_ages_all_samples,
                    ahe_sample_depths,
                    ahe_ages_all_samples_SE):

            #show AHe ages:
            depths = np.ones(len(ahe_ages_sample)) * ahe_sample_depth
            leg_data = ax_ahe.errorbar(ahe_ages_sample, depths,
                                       xerr=ahe_ages_sample_SE * 1.96,
                                       **erb_props)
        
        # show violin plots for AHe pdf
        for ahe_age_pdf, ahe_sample_depth in \
                zip(ahe_age_pdfs, ahe_sample_depths):

            # violin plots of single grain age pdf
            violin_width = max_depth / 20.0
            pdf_threshold = 1e-5

            ahe_age_pdf_combined = np.zeros_like(ahe_age_pdf[0])
            for pdf in ahe_age_pdf:
                ahe_age_pdf_combined += pdf

            ind_vp = ahe_age_pdf_combined > pdf_threshold

            vd = dict(coords=ahe_age_bin[ind_vp],
                      vals=ahe_age_pdf_combined[ind_vp],
                      mean=1.0, min=1.0, max=1.0, median=1.0)
            vp = ax_ahe.violin([vd],
                               positions=[ahe_sample_depth],
                               vert=False,
                               widths=violin_width,
                               showextrema=False)
            for pc in vp['bodies']:
                pc.set_edgecolor('darkblue')
                pc.set_facecolor('lightblue')

        data_label.append('He ages')

    # scatter plots of modelled vs observed ages
    if show_age_scatter and len(ax_aft_scatter_list) > 0:
        cmap_sc = pl.get_cmap('tab10')
        for rank, (j, ax_sc) in enumerate(zip(aft_sg_sorted, ax_aft_scatter_list)):
            sg_ages = np.array(single_grain_aft_ages[j])
            sg_se_min = np.array(single_grain_aft_ages_se_min[j])
            sg_se_plus = np.array(single_grain_aft_ages_se_plus[j])
            sim_min = float(aft_data_samples['simulated_AFT_min'].iloc[j])
            sim_max = float(aft_data_samples['simulated_AFT_max'].iloc[j])
            sim_mid = (sim_min + sim_max) / 2.0
            color = cmap_sc(rank % 10)
            y_lo = np.full(len(sg_ages), sim_mid - sim_min)
            y_hi = np.full(len(sg_ages), sim_max - sim_mid)
            ax_sc.errorbar(
                sg_ages, np.full(len(sg_ages), sim_mid),
                xerr=[sg_se_min * 1.96, sg_se_plus * 1.96],
                yerr=[y_lo, y_hi],
                fmt='o', ms=4, color=color, lw=0.75, capsize=2,
            )
            age_max_sc = max(float(sg_ages.max()), sim_max) * 1.05
            ax_sc.plot([0, age_max_sc], [0, age_max_sc], '--', color='grey', lw=0.8, zorder=0)
            ax_sc.set_xlim(0, age_max_sc)
            ax_sc.set_ylim(0, age_max_sc)
            ax_sc.set_box_aspect(1)
            ax_sc.set_title(f'{aft_sample_names[j]}', loc='left', fontsize="x-small")
            ax_sc.spines['top'].set_visible(False)
            ax_sc.spines['right'].set_visible(False)
        ax_aft_scatter_list[-1].set_xlabel('measured (Ma)')
        mid_idx = len(ax_aft_scatter_list) // 2
        ax_aft_scatter_list[mid_idx].set_ylabel('modelled (Ma)')

    if show_age_scatter and len(ax_ahe_scatter_list) > 0:
        cmap_sc = pl.get_cmap('tab10')
        for rank, (i, ax_sc) in enumerate(zip(ahe_sorted_indices, ax_ahe_scatter_list)):
            grain_ages = np.array(ahe_ages_all_samples[i])
            grain_ses = np.array(ahe_ages_all_samples_SE[i])
            mod_ages = np.array(modeled_ahe_age_samples[i])
            color = cmap_sc(rank % 10)
            all_vals = list(grain_ages)
            for g, g_age in enumerate(grain_ages):
                y_vals = mod_ages[g, :] if mod_ages.ndim == 2 else np.array([mod_ages[g]])
                y_center = float(np.mean(y_vals))
                y_lo = y_center - float(y_vals.min())
                y_hi = float(y_vals.max()) - y_center
                ax_sc.errorbar(
                    g_age, y_center,
                    xerr=grain_ses[g] * 1.96,
                    yerr=[[y_lo], [y_hi]],
                    fmt='o', ms=4, color=color, lw=0.75, capsize=2,
                )
                all_vals.extend(y_vals.tolist())
            age_max_sc = max(all_vals) * 1.05
            ax_sc.plot([0, age_max_sc], [0, age_max_sc], '--', color='grey', lw=0.8, zorder=0)
            ax_sc.set_xlim(0, age_max_sc)
            ax_sc.set_ylim(0, age_max_sc)
            ax_sc.set_box_aspect(1)
            depth = ahe_sample_depths[i]
            name = ahe_data_samples['sample'].values[i]
            ax_sc.set_title('%s (%.0f m)' % (name, depth), loc='left', fontsize="x-small")
            ax_sc.spines['top'].set_visible(False)
            ax_sc.spines['right'].set_visible(False)
        ax_ahe_scatter_list[-1].set_xlabel('measured (Ma)')
        mid_idx = len(ax_ahe_scatter_list) // 2
        ax_ahe_scatter_list[mid_idx].set_ylabel('modelled (Ma)')

    # add labels:
    axb.set_ylabel('Burial depth (m)')

    if contour_variable == 'salinity':
        axst.set_ylabel('Salinity\ntop bnd\n(kg/kg)')
        axhf.set_ylabel('Salinity\nlower bnd\n(kg/kg)')
    elif contour_variable == 'pressure':
        # no top or lower boundary line is plotted for pressure (see
        # above), so these panels are left without a y-axis label
        pass
    else:
        axst.set_ylabel('Surface\nT (%sC)' % degree_symbol)
        axhf.set_ylabel(r'HF (mW m$^{-2}$)', labelpad=12)

    axhf.set_xlabel('Time (Ma)')
    if show_temp_panel:
        ax_temp.set_xlabel('T (%sC)' % degree_symbol)

    if show_porosity_panel:
        ax_porosity.set_xlabel('Porosity (-)')

    if show_pressure_panel:
        ax_pressure.set_xlabel('Pressure (MPa)')

    if C_data is not None:
        ax_c.set_xlabel('Salinity (kg/kg)')
    if VR_model_data is not None:
        ax_vr.set_xlabel('VR (Ro)')
    if AFT_data is not None:
        ax_afta.set_xlabel('AFT age (Ma)')
    if AHe_data is not None:
        ax_ahe.set_xlabel('He age (Ma)')
    #ax_aftln.set_xlabel(r'AFT ln ($\mu m$)')

    if show_strat_column is True:
        ax_strat.set_xticks([])

    for ax in all_panels[3:]:
        ax.set_yticklabels([])

    for ax in [axst, axb]:
        ax.set_xticklabels([])

    for ax in all_panels:
        ax.yaxis.grid(False)
        ax.xaxis.grid(False)
        #ax.spines['right'].set_color('none')
        #ax.spines['top'].set_color('none')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

    #
    max_time = time_array_bp.max() / 1e6 * 1.1

    if (AFT_data is not None or AHe_data is not None or show_thermochron_data is False) \
            and show_prov_ages_simple is True:
        max_time = np.array(prov_ages).max() * 1.1

    if (AFT_data is not None and show_provenance_hist is True
            and simulated_AFT_data is not None):
        start_times = np.array([ai[0]
                                for a in aft_node_times_burial
                                for ai in a])
        max_time = start_times.max() * 1.1

    if (AHe_data is not None and simulated_AHe_data is not None and show_provenance_hist is True):
        start_times = np.array([ai[0]
                                for a in ahe_node_times_burial
                                for ai in a])
        max_time = start_times.max() * 1.1

    if max_age_burial_panel is not None:
        max_time_b = max_age_burial_panel
    else:
        max_time_b = max_time
    for ax in time_panels:
        ax.set_xlim(max_time_b, 0)

    for ax in depth_panels:
        ax.set_ylim(max_depth, -max_depth / 20.0)

    max_T = T_nodes[-1][active_nodes[-1]].max()

    if T_data is not None:
        max_T2 = T_obs.max()
        max_T = max(max_T, max_T2)
        
    if show_temp_panel:
        ax_temp.set_xlim(0, max_T * 1.2)

    if show_porosity_panel:
        max_porosity = porosity_nodes[-1][active_nodes[-1]].max()
        if porosity_data is not None and len(porosity_data['data']) > 0:
            max_porosity = max(max_porosity, porosity_data['data']['porosity'].max())
        ax_porosity.set_xlim(0, max_porosity * 1.2)

    if show_pressure_panel:
        max_pressure = total_pressure.max()
        if lithostatic_pressure is not None:
            max_pressure = max(max_pressure, lithostatic_pressure.max())
        if pressure_data is not None and len(pressure_data['data']) > 0:
            max_pressure = max(max_pressure, pressure_data['data']['FSIP_MPa'].max())
        ax_pressure.set_xlim(0, max_pressure * 1.1)

    if contour_variable == 'salinity':
        max_C = C_nodes[-1].max()

    if C_data is not None and len(salinity_data) > 0:
        if salinity_data.max() > max_C:
            max_C = salinity_data.max()
        ax_c.set_xlim(0, max_C * 1.1)

    if VR_model_data is not None:
        if vr_nodes is not None:
            max_VR = vr_nodes.max()
        else:
            max_VR = 1.5
        if len(vr_obs) > 0 and vr_obs.max() > max_VR:
            max_VR = vr_obs.max()
        ax_vr.set_xlim(0.1, max_VR * 1.1)

    thermochron_age_max = max_time

    if AFT_data is not None:
        if simulated_AFT_data is not None:
            thermochron_age_max = aft_age_nodes[active_nodes[-1]].max()
        if len(aft_age) > 0 and aft_age.max() > thermochron_age_max:
            thermochron_age_max = aft_age.max()

    if AHe_data is not None:
        if simulated_AHe_data is not None:
            if np.max(np.array(ahe_age_nodes)) > thermochron_age_max:
                thermochron_age_max = np.max(np.array(ahe_age_nodes))
        for ahe_ages_sample in ahe_ages_all_samples:
            if ahe_ages_sample.max() > thermochron_age_max:
                thermochron_age_max = ahe_ages_sample.max()
        if modeled_ahe_age_samples_max is not None:
            for mm in modeled_ahe_age_samples_max:
                if np.max(mm) > thermochron_age_max:
                    thermochron_age_max = np.max(mm)
                    logger.info('updated thermochron max age to %.2f from modeled AHe ages' % thermochron_age_max)
                else:
                    logger.info('modeled AHe ages max %.2f not changing thermochron max age %.2f' % (np.max(mm), thermochron_age_max))

    if max_age_thermochron_panel is not None:
        thermochron_age_max = max_age_thermochron_panel

    if AFT_data is not None:
        ax_afta.set_xlim(thermochron_age_max * 1.1, 0)
    if AHe_data is not None:
        ax_ahe.set_xlim(thermochron_age_max * 1.1, 0)

    # connecting lines from sample depth in age panel to scatter panel top/bottom
    if show_age_scatter and AFT_data is not None and len(ax_aft_scatter_list) > 0:
        x_right = ax_afta.get_xlim()[1]
        for j, ax_sc in zip(aft_sg_sorted, ax_aft_scatter_list):
            depth = aft_age_depth[j]
            for y_frac in (1.0, 0.0):
                con = mpatches.ConnectionPatch(
                    xyA=(x_right, depth), xyB=(0.0, y_frac),
                    coordsA='data', coordsB='axes fraction',
                    axesA=ax_afta, axesB=ax_sc,
                    color='lightgrey', lw=0.75, clip_on=False,
                )
                fig.add_artist(con)

    if show_age_scatter and AHe_data is not None and len(ax_ahe_scatter_list) > 0:
        x_right = ax_ahe.get_xlim()[1]
        for i, ax_sc in zip(ahe_sorted_indices, ax_ahe_scatter_list):
            depth = ahe_sample_depths[i]
            for y_frac in (1.0, 0.0):
                con = mpatches.ConnectionPatch(
                    xyA=(x_right, depth), xyB=(0.0, y_frac),
                    coordsA='data', coordsB='axes fraction',
                    axesA=ax_ahe, axesB=ax_sc,
                    color='lightgrey', lw=0.75, clip_on=False,
                )
                fig.add_artist(con)

    #ax_aftln.set_xlim(2, 17)
    #if max_T > 75.0:
    #    t_ticks = np.arange(0.0, max_T + 25.0, 25.0)
    #    ax_temp.set_xticks(t_ticks)

    # remove last tick label to avoid overlap
    if show_temp_panel:
        ax_temp.set_xticks(ax_temp.get_xticks()[:-1])
    if VR_model_data is not None:
        ax_vr.set_xticks(ax_vr.get_xticks()[:-1])
    if AFT_data is not None:
        ax_afta.set_xticks(ax_afta.get_xticks()[:-1])
    if AHe_data is not None:
        ax_ahe.set_xticks(ax_ahe.get_xticks()[:-1])

    for ax in all_panels[3:]:
        # reduce number of tick labels
        ax.set_xticks(ax.get_xticks()[::2])

    if contour_variable == 'salinity':
        axst.set_yticks(axst.get_yticks()[1::2])
        axhf.set_yticks(axhf.get_yticks()[::2])
    else:
        #hf_min = int(np.floor(basal_hf_array.min() * 100.0)) * 10.0
        #hf_max = int(np.ceil(basal_hf_array.max() * 100.0)) * 10.0
        #hf_ticks = np.arange(hf_min, hf_max + 5.0, 5.0)
        #axhf.set_yticks(hf_ticks)

        axhf.set_yticks(axhf.get_yticks()[::3])

        st_min = int(np.floor(surface_temp_array.min() / 10.0)) * 10.0
        st_max = int(np.ceil(surface_temp_array.max() / 10.0)) * 10.0
        st_ticks = np.arange(st_min, st_max + 5.0, 5.0)
        axst.set_yticks(st_ticks)

    if (show_temp_panel and T_data is not None and np.isnan(T_gof) == False
            and 'temperature' in show_gof_stats):
        ax_temp.text(0.5, 1.03,
                     'GOF=%0.2f\nRMSE=%0.1f' % (T_gof, T_rmse),
                     transform=ax_temp.transAxes,
                     **textprops)

    if (VR_model_data is not None and np.isnan(vr_GOF) == False
            and 'VR' in show_gof_stats):
        ax_vr.text(0.5, 1.03,
                   'GOF=%0.2f\nRMSE=%0.2f' % (vr_GOF, vr_rmse),
                   transform=ax_vr.transAxes,
                   **textprops)

    if (AFT_data is not None and np.isnan(aft_age_GOF) == False
            and 'AFT' in show_gof_stats):
        ax_afta.text(0.5, 1.03,
                     'GOF=%0.2f\nerror=%0.2f My'
                     % (aft_age_GOF, aft_age_error),
                     transform=ax_afta.transAxes,
                     **textprops)

    if (AHe_data is not None and np.isnan(ahe_age_gof) == False
            and 'AHe' in show_gof_stats):
        ax_ahe.text(0.5, 1.03,
                    'GOF=%0.2f\nerror=%0.2f My'
                    % (ahe_age_gof, ahe_age_error),
                    transform=ax_ahe.transAxes,
                    **textprops)

    if (show_pressure_panel and pressure_data is not None
            and np.isnan(pressure_gof) == False
            and 'pressure' in show_gof_stats):
        ax_pressure.text(0.5, 1.03,
                         'GOF=%0.2f\nRMSE=%0.1f' % (pressure_gof, pressure_rmse),
                         transform=ax_pressure.transAxes,
                         **textprops)

    #gs.tight_layout(fig, h_pad=0.02, w_pad=0.02)
    # add colorbar
    if len(depth_panels) > 1:
        cax_left = depth_panels[1].get_position().x0 + cb_buffer_hor
    else:
        cax_left = axb.get_position().x1 + cb_buffer_hor
    pos_right = all_panels[-1].get_position()
    cax_right = pos_right.x0 + pos_right.width
    cax_width = cax_right - cax_left - cb_buffer_hor
    cax_bottom = axhf.get_position().y0 + cb_buffer_vert
    cax = fig.add_axes([cax_left, cax_bottom, cax_width, 0.015])
    cb = fig.colorbar(tc, cax=cax, orientation='horizontal')
    cb.set_label(cb_label, fontsize='medium')

    for p in all_panels:
        locy = ticker.MaxNLocator(nbins=3)  # this locator puts ticks at regular intervals
        locx = ticker.MaxNLocator(nbins=3)  # this locator puts ticks at regular intervals

        p.xaxis.set_major_locator(locx)
        p.yaxis.set_major_locator(locy)

    if show_thermochron_data is False and contour_variable != 'salinity':
        # fewer ticks in colorbar in case of small space
        max_T = cb.locator().max()
        max_T_tick = np.ceil(max_T / 50.0) * 50.0
        T_ticks = np.arange(0, max_T_tick+50, 50.0)
        cb.set_ticks(T_ticks)

    tick_locator = ticker.MaxNLocator(nbins=4)
    cb.locator = tick_locator
    cb.update_ticks()

    model_label_merged = 'modeled ' + ', '.join(model_label)
    model_range_label_merged = 'modeled range ' + ', '.join(model_range_label)
    data_label_merged = 'observed ' + ', '.join(data_label)

    if add_legend is True:
        if leg_data is not None:
            leg_items += [leg_data]
            leg_labels += [data_label_merged]

        if len(leg_data_ext) >= 1:
            leg_items += leg_data_ext
            leg_labels += data_ext_label

        leg_items += [leg_model]
        leg_labels += [model_label_merged]

        if leg_model_range is not None:
            leg_items += [leg_model_range]
            leg_labels += [model_range_label_merged]

        if leg_ahe_sample_max is not None:
            leg_items += [leg_ahe_sample_min]
            leg_labels += ['modeled He age range samples']

        fig.legend(leg_items, leg_labels,
                   loc='outside lower center', ncol=ncols_legend, fontsize=legend_fontsize,
                   frameon=False, numpoints=1, handlelength=2)

    if add_panel_titles is True:
        if panel_title_numbers is True:
            panel_labels_init = ['1', '2', '3', '4', '5',
                                 '6', '7', '8', '9', '10', '11']

        else:
            panel_labels_init = ['a', 'b', 'c', 'd', 'e',
                                 'f', 'g', 'h', 'i', 'j', 'k']

        panel_labels = ['%s%s' % (panel_title_prefix, p) for p in panel_labels_init]

        for panel, label in zip(all_panels, panel_labels):
            panel.text(0.03, 1.02, label,
                       horizontalalignment='left',
                       verticalalignment='top',
                       weight='extra bold',
                       transform=panel.transAxes,
                       fontsize=panel_label_fs)

    return fig
