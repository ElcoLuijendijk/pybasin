__author__ = 'elcopone'

import pdb
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib.mlab
import matplotlib.patches as mpatches

#import useful_functions

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
        print 'figure exceeding b5 paper size'

    # initialize figure
    if landscape is True:
        print 'landscape figure'
        xs = height
        ys = width
    else:
        print 'portrait figure'
        xs = width
        ys = height

    if units != 'inch':
        xs = xs / 25.4
        ys = ys / 25.4

    print 'init fig,  size = %0.1f x %0.1f inch' %(xs, ys)

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

    if fontsize_legend != None:
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
                         contour_variable='temperature',
                         add_legend=True, debug=False):

    """
    create a figure comparing 1D burial and thermal model results
    with vitrinite reflectance, apatite fission track and present-day
    temperature data
    """

    [time_array_bp,
     surface_temp_array, basal_hf_array,
     z_nodes, active_nodes, T_nodes,
     node_strat, node_age,
     T_data, C_data, VR_model_data, AFT_data, AHe_data] = \
        model_run_data

    if T_data is not None:
        (T_depth,
         T_obs,
         T_obs_sigma,
         T_gof, T_rmse) = T_data

    if C_data is not None:
        [C_nodes, surface_salinity_array, salinity_lwr_bnd,
         salinity_depth, salinity_data, salinity_data_unc,
         salinity_RMSE, q_solute_bottom, q_solute_top] = C_data

    if VR_model_data is not None:
        [vr_nodes,
         vr_depth,
         vr_obs,
         vr_min,
         vr_max,
         vr_obs_sigma,
         vr_GOF] = VR_model_data

    if AFT_data != None:
        [simulated_AFT_data,
         aft_sample_names,
         aft_age_depth,
         aft_age,
         aft_age_stderr_min,
         aft_age_stderr_plus,
         aft_length_mean,
         aft_length_std,
         aft_data_type,
         aft_age_samples,
         single_grain_aft_ages,
         single_grain_aft_ages_se_min,
         single_grain_aft_ages_se_plus,
         aft_age_bins,
         aft_age_pdfs,
         aft_age_GOF,
         aft_age_error,
         aft_sample_times,
         aft_sample_temps] = AFT_data

    if AHe_data is not None:
        [ahe_sample_depths,
         ahe_ages_all_samples,
         ahe_ages_all_samples_SE,
         ahe_age_bin,
         ahe_age_pdfs,
         modeled_ahe_age_samples,
         modeled_ahe_age_samples_min,
         modeled_ahe_age_samples_max,
         ahe_age_gof, ahe_age_error,
         simulated_AHe_data] = AHe_data

    nt_total, n_nodes = T_nodes.shape

    if AFT_data is not None and simulated_AFT_data is not None:
        (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
         aft_ln_mean_nodes, aft_ln_std_nodes,
         aft_node_times_burial, aft_node_zs,
         aft_node_times, aft_node_temps) = simulated_AFT_data

        _, n_prov_scenarios, n_kinetic_scenarios = aft_age_nodes.shape

    if AHe_data is not None and simulated_AHe_data is not None:
        (ahe_age_nodes, ahe_age_nodes_min, ahe_age_nodes_max,
         ahe_node_times_burial, ahe_node_zs) = simulated_AHe_data

    degree_symbol = unichr(176)

    xsize = 170.0 / 25.4
    golden_ratio = (1.0 + np.sqrt(5))/2.0
    ysize= xsize / golden_ratio

    if add_legend is True:
        ysize = ysize * 1.1

    fig = pl.figure(figsize=(xsize, ysize))

    font = {'family': 'sans-serif',
            'size': 9}

    pl.rc('font', **font)

    width_ratios = [8, 3]

    nrows = 3
    ncols = 2

    bottom = 0.12
    if add_legend is True:
        bottom = 0.23

    leg_items = []
    leg_labels = []
    model_label = []
    model_range_label = []
    data_label = []
    data_ext_label = []
    leg_data = None
    leg_data_ext = []
    leg_model_range = None

    max_depth = z_nodes.max() * 1.1

    # skip VR, AFT and AHe panels if no data
    if VR_model_data is not None and len(vr_obs) == 0:
        VR_model_data = None
    if AFT_data is not None and len(aft_age) == 0:
        AFT_data = None
    if AHe_data is not None and len(ahe_ages_all_samples) == 0:
        AHe_data = None

    if C_data is not None:
        C_panel_ind = ncols
        ncols += 1
        width_ratios.append(3)

    if VR_model_data is not None:
        vr_panel_ind = ncols
        ncols += 1
        width_ratios.append(3)

    if AFT_data is not None:
        aft_panel_ind = ncols
        ncols += 1
        width_ratios.append(4)

    if AHe_data is not None:
        ahe_panel_ind = ncols
        ncols += 1
        width_ratios.append(4)

    gs = gridspec.GridSpec(nrows, ncols,
                           wspace=0.06, hspace=0.08,
                           bottom=bottom, top=0.96, left=0.12, right=0.97,
                           width_ratios=width_ratios,
                           height_ratios=[1, 4, 1])

    axb = fig.add_subplot(gs[1, 0])
    axst = fig.add_subplot(gs[0, 0])
    axhf = fig.add_subplot(gs[2, 0])
    #ax_strat = fig.add_subplot(gs[1, 1])
    ax_temp = fig.add_subplot(gs[1, 1])

    all_panels = [axb, axhf, axst, ax_temp]

    if C_data is not None:
        ax_c = fig.add_subplot(gs[1, C_panel_ind])
        all_panels.append(ax_c)
    if VR_model_data is not None:
        ax_vr = fig.add_subplot(gs[1, vr_panel_ind])
        all_panels.append(ax_vr)
    if AFT_data is not None:
        ax_afta = fig.add_subplot(gs[1, aft_panel_ind])
        all_panels.append(ax_afta)
    if AHe_data is not None:
        ax_ahe = fig.add_subplot(gs[1, ahe_panel_ind])
        all_panels.append(ax_ahe)

    #ax_aftln = fig.add_subplot(gs[1, 4])

    depth_panels = [all_panels[0]] + all_panels[3:]
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
                 "color": 'gray',
                 "mec": 'black',
                 "mfc": 'gray',
                 "lw": 0.5,
                 "zorder": 10}

    textprops = {"fontsize": 'small',
                 'ha': 'center',
                 'va': 'bottom',
                 'weight': 'normal',
                 'bbox': dict(facecolor="white",
                              ec='white',
                              alpha=0.7)}

    provenance_color = 'darkgray'
    cmap = matplotlib.cm.get_cmap('jet')

    if contour_variable == 'salinity':
        cnt_var = C_nodes
        cnt_step = 0.005
        cb_label = 'salinity (kg/kg)'
    else:
        cnt_var = T_nodes
        if T_nodes.max() < 50:
            cnt_step = 2.5
        elif T_nodes.max() < 100:
            cnt_step = 5.0
        else:
            cnt_step = 10.0
        cb_label = 'T (%s C)' % degree_symbol

    # plot surface temperature
    if contour_variable == 'salinity':
        axst.plot(time_array_bp / 1e6, surface_salinity_array,
                  **line_props)
    else:
        axst.plot(time_array_bp / 1e6, surface_temp_array,
                  **line_props)

    ts = 1.0e5

    if z_nodes.max() < 1000:
        ys = 1.0
    else:
        ys = 10.0

    yi = np.arange(z_nodes.min(), z_nodes.max() + ys, ys)

    cnt_var_mask = cnt_var.copy()
    cnt_var_mask[active_nodes == False] = np.nan

    time_2d = np.zeros([nt_total, n_nodes])

    for nn in xrange(n_nodes):
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
    ind = act == True

    print 'gridding T or salinity data vs time'
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
        try:
            zi[:, -tsi] = z_interpolated
        except IndexError:
            pdb.set_trace()

    if gridding_ok is True:
        # find max depth at each timestep
        max_depth_time = np.max(z_nodes, axis=1)
        max_depth_time2 = np.interp(xi,
                                    (time_array_bp/1.0e6)[::-1],
                                    max_depth_time[::-1])

        # filter interpolated values that are deeper than deepest fm.
        for nti in xrange(len(xi)):
            zi.mask[yi > max_depth_time2[nti], nti] = True

        print 'color mesh:'
        #tc = axb.pcolormesh(xg, yg, zi2, cmap='jet')
        c_int = np.arange(0.0, cnt_var.max()+cnt_step, cnt_step)
        tc = axb.contourf(xi, yi, zi, c_int, cmap='jet', zorder=1.0)

    else:
        plot_int = 1
        tc = axb.scatter(x[ind][::plot_int],
                         y[ind][::plot_int],
                         c=z[ind][::plot_int],
                         edgecolor="black", lw=0.1,
                         s=10,
                         cmap=cmap)

    major_strat = [n[:4] for n in node_strat]
    strat_transition = [m != n for m, n in zip(major_strat[:-1],
                                               major_strat[1:])]
    strat_transition.append(True)

    if (AFT_data is not None or AHe_data is not None) and show_provenance_hist is True:

        if AFT_data is not None:
            burial = aft_node_times_burial
            depths = aft_node_zs
        else:
            burial = ahe_node_times_burial
            depths = ahe_node_zs


        strat_count = 0
        for xb, yb, strat_trans in zip(burial, depths,
                                       strat_transition):
            if strat_trans is True:
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

                xf = np.concatenate((xb[min_prov_ind], xb[max_prov_ind][::-1]))
                yf = np.concatenate((yb[min_prov_ind], yb[max_prov_ind][::-1]))

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
        for i in xrange(n_strat_trans):
            leg_strat_unit, = axb.plot(time_array_bp / 1e6, z_nodes[:, ind][:, i],
                                  color='black', lw=0.5, zorder=100)

        leg_items += [leg_strat_unit]
        leg_labels += ['major stratigraphic unit']

    # plot basal heat flow
    if contour_variable == 'salinity':
        axhf.axhline(y=salinity_lwr_bnd, **line_props)
    else:
        axhf.plot(time_array_bp / 1e6, basal_hf_array * 1000.0,
                  **line_props)
        axhf.set_ylim(basal_hf_array.min() * 1000.0 * 0.95,
                      basal_hf_array.max() * 1000.0 * 1.05)

    # plot surface temperature
    leg_model, = ax_temp.plot(T_nodes[-1, active_nodes[-1]],
                              z_nodes[-1, active_nodes[-1]],
                              **line_props)
    model_label.append('temperature')

    if T_data is not None and len(T_data) > 0:
        leg_data = ax_temp.errorbar(T_obs, T_depth,
                                    xerr=T_obs_sigma * 2, **erb_props)
        data_label.append('temperature')

    # plot modeled salinity
    if C_data is not None and C_nodes is not None:
        leg_model, = ax_c.plot(C_nodes[-1, active_nodes[-1]],
                               z_nodes[-1, active_nodes[-1]],
                               **line_props)
        model_label.append('salinity')

    #pdb.set_trace()
    if C_data is not None and len(salinity_data) > 0:
        leg_data = ax_c.scatter(salinity_data, salinity_depth,
                                **scatter_props)
        data_label.append('salinity')

    # plot vitrinite
    if VR_model_data is not None and vr_nodes is not None and len(vr_obs) > 0:
        leg_model, = ax_vr.plot(vr_nodes[-1, active_nodes[-1]],
                                z_nodes[-1, active_nodes[-1]],
                                **line_props)

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
        violin_width = max_depth / 20.0
        pdf_threshold = 1e-5
        for sample_no in xrange(len(aft_age)):
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
                    pc.set_edgecolor('darkblue')
                    pc.set_facecolor('lightblue')

        leg_violin = mpatches.Patch(facecolor='lightblue',
                                    edgecolor='blue')
        leg_data_ext.append(leg_violin)
        data_ext_label.append('age distribution')

        # show single grain AFT ages, without errorbar
        for sample_no in xrange(len(single_grain_aft_ages)):
            x = single_grain_aft_ages[sample_no]

            if x is not None:
                y = np.ones_like(x) * aft_age_depth[sample_no]
                leg_sg = ax_afta.scatter(x, y, color='black', s=5, marker='o')
                if len(leg_labels) == 0 or 'single grain' not in leg_labels[-1]:
                    leg_data_ext.append(leg_sg)
                    data_ext_label.append('single grain AFT ages')

        ind_ca = np.array([a is None for a in single_grain_aft_ages])

        if True in ind_ca:
            # show central ages
            leg_data = ax_afta.errorbar(aft_age[ind_ca], aft_age_depth[ind_ca],
                                        xerr=[aft_age_stderr_min[ind_ca] * 1.96,
                                              aft_age_stderr_plus[ind_ca] * 1.96],
                                        **erb_props)
            #if len(leg_labels) == 0 or 'AFT age' not in leg_labels[-1]:
            data_label.append('AFT age')

            #leg_items.append(leg_ca)
            #leg_labels.append('single grain AFT ages')

    if AFT_data is not None and simulated_AFT_data is not None:
        for n_prov in xrange(n_prov_scenarios):
            for n_kin in xrange(n_kinetic_scenarios):
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

        (ahe_age_nodes, ahe_age_nodes_min, ahe_age_nodes_max,
                     ahe_node_times_burial, ahe_node_zs) = simulated_AHe_data
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
        for n_prov in xrange(n_prov_scenarios):
            for n_rad in xrange(n_grain_radius):
                leg_model, = ax_ahe.plot(ahe_age_nodes_array[active_nodes[-1], n_rad, n_prov],
                                         z_nodes[-1, active_nodes[-1]],
                                         **line_props)

        model_label.append('AHe ages')
        model_range_label.append('AHe ages')


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

        data_label.append('AHe ages')

    # add labels:
    axb.set_ylabel('Burial depth (m)')

    if contour_variable == 'salinity':
        axst.set_ylabel('Salinity\ntop bnd\n(kg/kg)')
        axhf.set_ylabel('Salinity\nlower bnd\n(kg/kg)')
    else:
        axst.set_ylabel('Surface\nT (%sC)' % degree_symbol)
        axhf.set_ylabel(r'HF (mW m$^{-2}$)')

    axhf.set_xlabel('Time (Ma)')
    ax_temp.set_xlabel('T (%sC)' % degree_symbol)

    if C_data is not None:
        ax_c.set_xlabel('Salinity (kg/kg)')
    if VR_model_data is not None:
        ax_vr.set_xlabel('VR (Ro)')
    if AFT_data is not None:
        ax_afta.set_xlabel('AFT age (Ma)')
    if AHe_data is not None:
        ax_ahe.set_xlabel('AHe age (Ma)')
    #ax_aftln.set_xlabel(r'AFT ln ($\mu m$)')

    for ax in all_panels[3:]:
        ax.set_yticklabels([])

    for ax in [axst, axb]:
        ax.set_xticklabels([])

    for ax in all_panels:
        ax.yaxis.grid(True)
        ax.xaxis.grid(False)
        #ax.spines['right'].set_color('none')
        #ax.spines['top'].set_color('none')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

    #
    max_time = time_array_bp.max() / 1e6 * 1.1

    if (AFT_data is not None and show_provenance_hist is True
            and simulated_AFT_data is not None):
        start_times = np.array([ai[0]
                                for a in aft_node_times_burial
                                for ai in a])
        max_time = start_times.max() * 1.1

    for ax in time_panels:
        ax.set_xlim(max_time, 0)

    for ax in depth_panels:
        ax.set_ylim(max_depth, -max_depth / 20.0)

    max_T = T_nodes[-1].max()

    if T_data is not None:
        max_T = T_obs.max()

    ax_temp.set_xlim(0, max_T * 1.2)

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

    if AFT_data is not None:
        ax_afta.set_xlim(thermochron_age_max * 1.1, 0)
    if AHe_data is not None:
        ax_ahe.set_xlim(thermochron_age_max * 1.1, 0)

    #ax_aftln.set_xlim(2, 17)
    #if max_T > 75.0:
    #    t_ticks = np.arange(0.0, max_T + 25.0, 25.0)
    #    ax_temp.set_xticks(t_ticks)

    # remove last tick label to avoid overlap
    ax_temp.set_xticks(ax_temp.get_xticks()[:-1])
    if VR_model_data is not None:
        ax_vr.set_xticks(ax_vr.get_xticks()[:-1])
    if AFT_data is not None:
        ax_afta.set_xticks(ax_afta.get_xticks()[:-1])
    if AHe_data is not None:
        ax_ahe.set_xticks(ax_ahe.get_xticks()[:-1])

    for ax in all_panels[3:]:
        # reduce number of tick labels
        print ax.get_xticks()
        ax.set_xticks(ax.get_xticks()[::2])

    if contour_variable == 'salinity':
        axst.set_yticks(axst.get_yticks()[1::2])
        axhf.set_yticks(axhf.get_yticks()[::2])
    else:
        hf_min = int(np.floor(basal_hf_array.min() * 100.0)) * 10.0
        hf_max = int(np.ceil(basal_hf_array.max() * 100.0)) * 10.0
        hf_ticks = np.arange(hf_min, hf_max + 5.0, 5.0)
        axhf.set_yticks(hf_ticks)

        st_min = int(np.floor(surface_temp_array.min() / 10.0)) * 10.0
        st_max = int(np.ceil(surface_temp_array.max() / 10.0)) * 10.0
        st_ticks = np.arange(st_min, st_max + 5.0, 5.0)
        axst.set_yticks(st_ticks)

    if T_data is not None and np.isnan(T_gof) == False:
        ax_temp.text(0.5, 1.03,
                     'GOF=%0.2f' % T_gof,
                     transform=ax_temp.transAxes,
                     **textprops)

    if VR_model_data is not None and np.isnan(vr_GOF) == False:
        ax_vr.text(0.5, 1.03,
                   'GOF=%0.2f' % vr_GOF,
                   transform=ax_vr.transAxes,
                   **textprops)

    if AFT_data is not None and np.isnan(aft_age_GOF) == False:
        ax_afta.text(0.5, 1.03,
                     'GOF=%0.2f\nerror=%0.2f My'
                     % (aft_age_GOF, aft_age_error),
                     transform=ax_afta.transAxes,
                     **textprops)

    if AHe_data is not None and np.isnan(ahe_age_gof) == False:
        ax_ahe.text(0.5, 1.03,
                    'GOF=%0.2f\nerror=%0.2f My'
                    % (ahe_age_gof, ahe_age_error),
                    transform=ax_ahe.transAxes,
                    **textprops)

    #gs.tight_layout(fig, h_pad=0.02, w_pad=0.02)
    # add colorbar
    cax_left = ax_temp.get_position().x0 + 0.05
    pos_right = all_panels[-1].get_position()
    cax_right = pos_right.x0 + pos_right.width
    cax_width = cax_right - cax_left - 0.03
    cax_bottom = axhf.get_position().y0
    cax = fig.add_axes([cax_left, cax_bottom, cax_width, 0.015])
    cb = fig.colorbar(tc, cax=cax, orientation='horizontal')
    cb.set_label(cb_label, fontsize='medium')

    if contour_variable is 'salinity':
        cb_ticks = [0.0, 0.1, 0.2, 0.3, 0.4]
        cb.set_ticks(cb_ticks)

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

        fig.legend(leg_items, leg_labels,
                   loc='lower center', ncol=3, fontsize='small',
                   handlelength=3, frameon=False)

    return fig