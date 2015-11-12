__author__ = 'elcopone'

import pdb
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import matplotlib.mlab

#import useful_functions

def setup_figure( width = 125.0, height = 'g', fontsize = 'x-small',
                  landscape = False,  fontsize_legend = 0, units= 'mm'):

    '''
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

    '''

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
    if landscape  ==  True:
        print 'landscape figure'
        xs = height ; ys = width
    else:
        print 'portrait figure'
        xs = width ; ys = height

    if units != 'inch':
        xs = xs/25.4
        ys = ys/25.4

    print 'init fig,  size  =  %0.1f x %0.1f inch' %(xs, ys)

    fig  =  pl.figure(figsize = (xs,ys))

    # set default parameters for figure
    if type(fontsize)  ==  str:
        if fontsize  ==  'xxx-small':
            fontsize_s  =  'xx-small'
            fontsize_l  =  'xx-small'
            fontsize_leg  =  'xx-small'
        elif fontsize  ==  'xx-small':
            fontsize_s  =  'xx-small'
            fontsize_l  =  'x-small'
            fontsize_leg  =  'xx-small'
        elif fontsize  ==  'x-small':
            fontsize_s  =  'x-small'
            fontsize_l  =  'small'
            fontsize_leg  =  'xx-small'
        elif fontsize  ==  'small':
            fontsize_s  =  'small'
            fontsize_l  =  'medium'
            fontsize_leg  =  'x-small'
        elif fontsize  ==  'medium':
            fontsize_s  =  'medium'
            fontsize_l  =  'large'
            fontsize_leg  =  'small'
    else:
        fontsize_s  =  fontsize
        fontsize_l  =  'xx-small'
        fontsize_leg  =  fontsize

    if fontsize_legend != None:
        fontsize_leg  =  fontsize_legend

    params  =  {'axes.labelsize': fontsize_s,
                'text.fontsize': fontsize_l,
                'legend.fontsize': fontsize_leg,
                'axes.titlesize' : fontsize_l,
                'xtick.labelsize': fontsize_s,
                'ytick.labelsize': fontsize_s}

    pl.rcParams.update(params)

    return fig


def model_vs_data_figure(time_array_bp,
                         surface_temp_array, basal_hf_array,
                         z_nodes, active_nodes, T_nodes,
                         node_strat, node_age,
                         simulated_AFT_data,
                         vr_nodes,
                         T_depth, T_data, T_data_sigma,
                         vr_depth, vr_data, vr_data_sigma,
                         aft_age_depth, aft_age, aft_age_min_ci, aft_age_plus_ci,
                         aft_length_mean, aft_length_std, aft_data_type,
                         T_GOF, vr_GOF, age_GOF,
                         show_provenance_hist=True,
                         time_int_grid=10,
                         model_data_fig_bw=False):

    """
    create a figure comparing 1D burial and thermal model results
    with vitrinite reflectance, apatite fission track and present-day
    temperature data

    :param time_array_bp:
    :param surface_temp_array:
    :param basal_hf_array:
    :param z_nodes:
    :param active_nodes:
    :param T_nodes:
    :param node_strat:
    :param node_age:
    :param aft_node_times_burial:
    :param aft_node_zs:
    :param aft_age_nodes:
    :param aft_age_nodes_min:
    :param aft_age_nodes_max:
    :param aft_ln_mean_nodes:
    :param aft_ln_std_nodes:
    :param vr_nodes:
    :param T_depth:
    :param T_data:
    :param T_data_sigma:
    :param vr_depth:
    :param vr_data:
    :param vr_data_sigma:
    :param aft_age_depth:
    :param aft_age:
    :param aft_age_min_ci:
    :param aft_age_plus_ci:
    :param aft_length_mean:
    :param aft_length_std:
    :param aft_data_type:
    :param T_GOF:
    :param vr_GOF:
    :param age_GOF:
    :param show_provenance_hist:
    :param time_int_grid:
    :return:

    """

    nt_total, n_nodes = T_nodes.shape

    if simulated_AFT_data is not None:
        (aft_age_nodes, aft_age_nodes_min, aft_age_nodes_max,
         aft_ln_mean_nodes, aft_ln_std_nodes,
         aft_node_times_burial, aft_node_zs) = simulated_AFT_data

        _, n_prov_scenarios, n_kinetic_scenarios = aft_age_nodes.shape

    #n_nodes = T_nodes.shape[1]
    degree_symbol = unichr(176)

    fig = setup_figure(width='2col', height=0.5, fontsize='xx-small')

    gs = gridspec.GridSpec(3, 4,
                           wspace=0.0, hspace=0.0,
                           width_ratios=[8, 2, 2, 4],
                           height_ratios=[1, 5, 1])

    axb = fig.add_subplot(gs[1, 0])
    axst = fig.add_subplot(gs[0, 0])
    axhf = fig.add_subplot(gs[2, 0])
    #ax_strat = fig.add_subplot(gs[1, 1])
    ax_temp = fig.add_subplot(gs[1, 1])
    ax_vr = fig.add_subplot(gs[1, 2])
    ax_afta = fig.add_subplot(gs[1, 3])
    #ax_aftln = fig.add_subplot(gs[1, 4])

    all_panels = [axb, axst, axhf, ax_temp, ax_vr, ax_afta]
    depth_panels = [axb, ax_temp, ax_vr, ax_afta]
    time_panels = [axb, axhf, axst]

    line_props = {"color": "black", "lw": 1.0}
    sc_props = {"facecolor": "gray",
                "edgecolor": "black",
                "s": 30}
    erb_props = {"marker": "o",
                 "ms": 4,
                 "linestyle": "None",
                 "color": 'gray',
                 "mec": 'black',
                 "mfc": 'gray',
                 "zorder": 10}

    erb_props2 = {"marker": "^",
                  "ms": 3,
                  "linestyle": "None",
                  "color": 'gray',
                  "mec": 'gray',
                  "mfc": 'gray',
                  "zorder": 10}

    textprops = {"fontsize": 'xx-small',
                 'ha': 'center',
                 'va': 'bottom',
                 'weight': 'bold',
                 'bbox': dict(facecolor="white", ec='white', alpha=0.7)}

    if model_data_fig_bw is True:
        provenance_color = 'black'
        cmap = matplotlib.cm.get_cmap('Greys')

    else:
        provenance_color = 'darkgray'
        cmap = matplotlib.cm.get_cmap('jet')


    # plot surface temperature
    axst.plot(time_array_bp / 1e6, surface_temp_array,
              **line_props)

    # burial history and temperature
    time_int_grid1 = 50

    xi = (time_array_bp / 1e6)[::time_int_grid1]
    yi = np.linspace(z_nodes.min(), z_nodes.max(),
                     (z_nodes.max() - z_nodes.min()) / 10.0)

    T_nodes_mask = T_nodes.copy()
    T_nodes_mask[active_nodes == False] = np.nan

    time_2d = np.zeros([nt_total, n_nodes])

    for nn in xrange(n_nodes):
        time_2d[:, nn] = time_array_bp / 1.0e6

    x = time_2d[::time_int_grid].ravel()
    y = z_nodes[::time_int_grid].ravel()
    z = T_nodes_mask[::time_int_grid].ravel()
    act = active_nodes[::time_int_grid].ravel()

    ind = act == True

    plot_int = 1

    print 'gridding T vs time data'
    #zi = matplotlib.mlab.griddata(x[ind][::plot_int],
    #                              y[ind][::plot_int],
    #                              z[ind][::plot_int],
    #                              xi, yi)

    # find max depth
    max_depth_time = np.max(z_nodes, axis=1)

    # filter interpolated values that are deeper than deepest fm.
    #for nt in xrange(len(xi)):
    #    zi.mask[yi > max_depth_time[nt], nt] = True

    print 'color mesh:'
    #tc = axb.pcolormesh(xi, yi, zi, cmap='jet')
    c_int = np.arange(T_nodes.min(), T_nodes.max()+5.0, 5.0)

    #tc = axb.contourf(xi, yi, zi, c_int, cmap='jet')
    plot_int = 1
    tc = axb.scatter(x[ind][::plot_int],
                     y[ind][::plot_int],
                     c=z[ind][::plot_int],
                     edgecolor="None",
                     s=3,
                     cmap=cmap)

    major_strat = [n[:4] for n in node_strat]
    strat_transition = [m != n for m, n in zip(major_strat[:-1],
                                               major_strat[1:])]
    strat_transition.append(True)

    #figb = pl.figure()
    #axbt = figb.add_subplot(1, 1, 1)

    if show_provenance_hist is True and simulated_AFT_data is not None:
        for xb, yb, strat_trans in zip(aft_node_times_burial,
                                       aft_node_zs,
                                       strat_transition):
            if strat_trans is True:
                xf = np.concatenate((xb[0], xb[-1][::-1]))
                yf = np.concatenate((yb[0], yb[-1][::-1]))

                axb.fill(xf, yf, color='lightgrey', zorder=0)

                for xbi, ybi in zip(xb, yb):
                    axb.plot(xbi, ybi, color=provenance_color, lw=0.5)

    else:
        ind = np.array(strat_transition) == True
        axb.plot(time_array_bp / 1e6, z_nodes[:, ind], color='black', lw=0.5)

    # plot basal heat flow
    axhf.plot(time_array_bp / 1e6, basal_hf_array * 1000.0,
              **line_props)

    # plot temperature
    ax_temp.plot(T_nodes[-1, active_nodes[-1]],
                 z_nodes[-1, active_nodes[-1]],
                 **line_props)

    ax_temp.errorbar(T_data, T_depth, xerr=T_data_sigma * 2, **erb_props)

    # plot vitrinite
    if vr_nodes is not None:
        ax_vr.plot(vr_nodes[-1, active_nodes[-1]],
                   z_nodes[-1, active_nodes[-1]],
                   **line_props)

    ax_vr.errorbar(vr_data, vr_depth, xerr=vr_data_sigma, **erb_props)

    # plot aft ages
    if simulated_AFT_data is not None:
        ax_afta.fill_betweenx(z_nodes[-1, active_nodes[-1]],
                              aft_age_nodes_min[active_nodes[-1]],
                              aft_age_nodes_max[active_nodes[-1]],
                              color='lightgrey')

    ax_afta.plot(node_age[active_nodes[-1]], z_nodes[-1, active_nodes[-1]],
                 color='blue', lw=1.5, ls='---', zorder=101)

    ind_ca = aft_data_type == 'central_age'
    ind_pop = aft_data_type == 'age_population'

    #pdb.set_trace()

    #ax_afta.errorbar(aft_age[ind_pop], aft_age_depth[ind_pop],
    #                 xerr=[aft_age_plus_ci[ind_pop], aft_age_min_ci[ind_pop]],
    #                 **erb_props2)

    ax_afta.errorbar(aft_age[ind_ca], aft_age_depth[ind_ca],
                     xerr=[aft_age_plus_ci[ind_ca], aft_age_min_ci[ind_ca]],
                     **erb_props)

    if simulated_AFT_data is not None:
        for n_prov in xrange(n_prov_scenarios):
            for n_kin in xrange(n_kinetic_scenarios):
                ax_afta.plot(aft_age_nodes[active_nodes[-1], n_prov, n_kin],
                             z_nodes[-1, active_nodes[-1]],
                             **line_props)

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

    # add labels:
    axst.set_ylabel('Surface T (%sC)' % degree_symbol)
    axb.set_ylabel('Burial depth (m)')
    axhf.set_ylabel(r'HF (mW m$^{-2}$)')

    axhf.set_xlabel('Time (Ma)')
    ax_temp.set_xlabel('T (%sC)' % degree_symbol)
    ax_vr.set_xlabel('VR (Ro)')
    ax_afta.set_xlabel('AFT age (Ma)')
    #ax_aftln.set_xlabel(r'AFT ln ($\mu m$)')

    for ax in [ax_temp, ax_vr, ax_afta]:
        ax.set_yticklabels([])

    for ax in [axst, axb]:
        ax.set_xticklabels([])

    for ax in [axst, axb, axhf, ax_temp, ax_vr, ax_afta]:
        ax.grid()

    #
    max_depth = z_nodes.max() * 1.1
    max_time = time_array_bp.max() / 1e6 * 1.1

    if show_provenance_hist is True and simulated_AFT_data is True:
        start_times = np.array([ai[0]
                                for a in aft_node_times_burial
                                for ai in a])
        max_time = start_times.max() * 1.1

    for ax in time_panels:
        ax.set_xlim(max_time, 0)

    for ax in depth_panels:
        ax.set_ylim(max_depth, -20.0)

    max_T = T_nodes[-1].max()
    if vr_nodes is not None:
        max_VR = vr_nodes.max()
    else:
        max_VR = 1.5
    if simulated_AFT_data is not None:
        afta_max = aft_age_nodes[active_nodes[-1]].max()
    else:
        afta_max = max_time

    # TODO: this fails when no T data, fix this
    try:
        if T_data.max() > max_T:
            max_T = T_data.max()
    except ValueError:
        print 'no T data, continuing'

    try:
        if vr_data.max() > max_VR:
            max_VR = vr_data.max()
    except ValueError:
        print 'no VR data, continuing'

    if aft_age[ind_ca].max() > afta_max:
        afta_max = aft_age[ind_ca].max()

    ax_temp.set_xlim(0, max_T * 1.1)
    ax_vr.set_xlim(0.1, max_VR * 1.1)
    ax_afta.set_xlim(afta_max * 1.1, 0)
    #ax_aftln.set_xlim(2, 17)

    print ax_temp.get_xticks()

    for ax in [ax_vr, ax_afta]:
        if len(ax.get_xticks()[1::2]) >= 2:
            ax.set_xticks(ax.get_xticks()[1::2])
        else:
            ax.set_xticks(ax.get_xticks()[::2])

    ax_vr.set_xticks(ax_vr.get_xticks()[:-1])

    t_ticks = np.arange(25, ax_temp.get_xticks()[-1], 50.0)
    ax_temp.set_xticks(t_ticks)

    hf_min = int(np.floor(basal_hf_array.min() * 100.0)) * 10.0
    hf_max = int(np.ceil(basal_hf_array.max() * 100.0)) * 10.0
    hf_ticks = np.arange(hf_min, hf_max + 5.0, 5.0)
    axhf.set_yticks(hf_ticks)

    st_min = int(np.floor(surface_temp_array.min() / 10.0)) * 10.0
    st_max = int(np.ceil(surface_temp_array.max() / 10.0)) * 10.0
    st_ticks = np.arange(st_min, st_max + 5.0, 5.0)
    axst.set_yticks(st_ticks)

    print 'new ticks'
    print ax_temp.get_xticks()

    # add colorbar
    #cax = useful_functions.add_subplot_axes(axb, [0.01, 0.14, 0.5, 0.025])
    cax = fig.add_axes([0.65, 0.1, 0.25, 0.015])
    cb = fig.colorbar(tc, cax=cax, orientation='horizontal')
    cb.set_label('T (%s C)' % degree_symbol, fontsize='xx-small')

    if np.isnan(T_GOF) == False:
        ax_temp.text(0.5, 1.03,
                     'GOF=%0.2f' % T_GOF,
                     transform=ax_temp.transAxes,
                     **textprops)

    if np.isnan(vr_GOF) == False:
        ax_vr.text(0.5, 1.03,
                   'GOF=%0.2f' % vr_GOF,
                   transform=ax_vr.transAxes,
                   **textprops)

    ax_afta.text(0.5, 1.03,
                 'GOF=%0.2f' % age_GOF,
                 transform=ax_afta.transAxes,
                 **textprops)

    fig.tight_layout()

    return fig