try:
    from mpl_toolkits.basemap import Basemap
except:
    print 'no basemap present'
import numpy as np
#import matplotlib
#matplotlib.use('GTKagg')

import matplotlib.pyplot as pl
import sys,  math,  random,  os,  shutil,  pdb, string
import tables
import matplotlib.patches as patches
import matplotlib.mlab

sys.path.append("../Libraries/")
import plotLib


def initFigure(vert_size):
    # initialize figure
    pl.figure(figsize = (190./25.4, vert_size/25.4))  # set figure size to a4 
    # set default parameters for figure
    params  = {'axes.labelsize': 'xx-small', 'text.fontsize': 'xx-small'\
    , 'legend.fontsize': 'xx-small', 'text.fontsize': 'xx-small'\
    , 'xtick.labelsize': 'xx-small', 'axes.titlesize' : 'xx-small'\
    , 'ytick.labelsize': 'xx-small'}
    pl.rcParams.update(params)
    return

    
def initLandscapeFigure():
    
    golden_ratio = (1.0 + math.sqrt(5))/2.0
    width = 190.0
    length = width * golden_ratio
    
    # initialize figure
    pl.figure(figsize = (length/25.4, width/25.4))  # set figure size to a4 
    # set default parameters for figure
    params  = {'axes.labelsize': 'xx-small', 'text.fontsize': 'xx-small'\
    , 'legend.fontsize': 7, 'text.fontsize': 'xx-small'\
    , 'xtick.labelsize': 'xx-small', 'axes.titlesize' : 'xx-small'\
    , 'ytick.labelsize': 'xx-small'}
    pl.rcParams.update(params)
    return


def initLandscapeFigure_poster():
    
    golden_ratio = (1.0 + math.sqrt(5))/2.0
    width = 190.0
    length = width*golden_ratio
    
    # initialize figure
    pl.figure(figsize = (length/25.4, width/25.4))  # set figure size to a4 
    # set default parameters for figure
    params  = {'axes.labelsize': 'x-small', 'text.fontsize': 'xx-small'\
    , 'legend.fontsize': 'xx-small'\
    , 'xtick.labelsize': 'x-small', 'axes.titlesize' : 'xx-small'\
    , 'ytick.labelsize': 'x-small'}
    pl.rcParams.update(params)
    return


def find_AFT_locs(h5file, AFTwellNames):
    
    coordinateTable = h5file.root.wells.coordinates 
    #AFTwellNames = ['ALM-01', 'AND - 06', 'BKZ-01', 'BRAK-01', 'HSW-01', 'KDK-01', 'LOZ-01', 'NDW-01', 'SMG-01', 'SPC-01', 'VEH-01', 'WWK-01', 'WWN-01']    
    #AFTwellNames = ['BKZ-01', 'HSW-01', 'KDK-01', 'NDW-01', 'SMG-01', 'SPC-01', 'WWK-01']    


    x_AFT = np.zeros((len(AFTwellNames)))
    y_AFT = np.zeros((len(AFTwellNames)))   


    for i, well in zip(range(len(AFTwellNames)), AFTwellNames):
        
        condition = "short_name  == '%s'" %well
        x_AFT[i] = float([row['long'] for row in coordinateTable.where(condition)][0])
        y_AFT[i] = float([row['lat'] for row in coordinateTable.where(condition)][0])
    
    return x_AFT, y_AFT


def plotBurialHistory(strat_steps, age_base_all, cellSize_steps):
    
    #######################################
    # plot burial history
    #######################################
    
    stratMarkers, strat_ages, stratColors = setStratCodeMarkersandColors()
    
    stratMarkers = stratMarkers[:: - 1]
    stratColors = stratColors[:: - 1]
    
    Nsteps = len(strat_steps)
    stratPlotDepths = np.zeros((len(stratMarkers), Nsteps)) ; stratPlotAges = np.zeros((Nsteps))
    
    # go through all tectonic time steps
    for i in xrange(Nsteps):
        stratPlotAges[i] = age_base_all[i]
        
        
        # sample the stratigraphic units present at each time step
        for j, stratName in enumerate(strat_steps[i]):
                                            
            # find the stratigraphic units in the stratigraphic marker list
            for h, stratIndex in enumerate(stratMarkers):
                a = len(stratIndex)
                
                # if stratigrpahic unit is present in marker list: store depth of base of stratigraphic unit
                if stratName[:a] == stratIndex:
                    #print 'plot strat %s' %stratIndex
                    stratPlotDepths[h, i] = cellSize_steps[i, :j].sum()
                elif type(stratIndex) == list:
                    #print 'plot compiled strat:'
                    for stratIndex_item in stratIndex:
                        a = len(stratIndex_item)
                            
                        if stratName[:a] == stratIndex_item:
                            #print stratIndex_item
                            stratPlotDepths[h, i] = cellSize_steps[i][:j].sum()
    
    # plot subsidence:
    for h in xrange(len(stratMarkers) - 1):
        pl.fill_between(stratPlotAges, (stratPlotDepths[h, :]/1000.), (stratPlotDepths[h + 1, :]/1000.), color = stratColors[h], label = stratMarkers[h])
        #pl.plot(stratPlotAges, (stratPlotDepths[h, :]/1000.))
    pl.fill_between(stratPlotAges, stratPlotDepths[-1, :]/1000., 0, color = stratColors[-1], label = stratMarkers[-1])
    
    return


def calculateDepths_mid(cellSize):
    
    # calculate cell depth midpoints from an 1D array of cell dimensions
    try:
        Ntectonic_steps, Ncells = cellSize.shape
        transformArray = False
    except:
        # transform 1D array to 2D:
        cellSize.shape = 1,  - 1
        Ntectonic_steps, Ncells = cellSize.shape
        transformArray = True
    
    cellDepths_mid = np.zeros((Ntectonic_steps, Ncells))
    cellDepths_mid[:, 0] = cellSize[:, 0]/2.
    for i in xrange(Ntectonic_steps):
        for j in range(1, Ncells):
                #cellDepths[j] = cellDepths[j - 1] + (cellSize[j - 1]/2 + cellSize[j]/2)
                cellDepths_mid[i, j] = cellSize[i, :j].sum() + (cellSize[i, j]/2.)
    
    if transformArray:
        # transform arrays back to 1D
        #print cellDepths_mid
        cellSize.shape = - 1
        cellDepths_mid.shape = - 1
    return cellDepths_mid


def plotTemperatureGrid(cellSize, age_base_all, temperature):
    
    ########################
    # plot temperautre grid:
    ########################
    #np.savetxt('temperatures.csv', temperature)
    
    # get array dimension
    xr, yr = temperature.shape
    
    #y = np.ravel(calculateDepths_mid(cellSize)/1000.0)
    #z = np.ravel(temperature)
    
    # create 2d array for timesteps:
    xg = np.zeros(temperature.shape)
    for i in xrange(xg.shape[0]):
        xg[i] = i
    
    # add extra row in arrays for surface temperatures:
    zp = np.concatenate((np.ones((xr, 1)), temperature), axis = 1)
    xp = np.concatenate((np.ones((xr, 1)), xg), axis = 1)
    yp = np.concatenate((np.zeros((xr, 1)), calculateDepths_mid(cellSize)/1000.0), axis = 1)
    zp[:, 0] = zp[:, 1]
    xp[:, 0] = xp[:, 1]
    
    # convert arrays to 1D column:
    x = np.ravel(xp) ; y = np.ravel(yp) ; z = np.ravel(zp)

    # create x and y coordinate arrays:
    ystep = 0.001
    #print 'y interval: ', ystep
    xi = np.arange(len(age_base_all))
    yi = np.arange(0, y.max() + ystep, ystep)

    # create array to deliniate grid
    #y_formask = np.ravel(cellSize)
    #ygrid = matplotlib.mlab.griddata(x, y, y_formask, xi, yi)
    zi = matplotlib.mlab.griddata(x, y, z, xi, yi)
    zi_masked = np.ma.masked_where(zi < 10, zi)
    
    V = np.arange(0, 120, 10)
    
    cft = pl.contourf(age_base_all, yi, zi_masked, V, cmap = pl.cm.jet)
    
    #cft = pl.imshow(zi_masked, extent = (xi.min(), xi.max(), yi.min(), yi.max()), aspect = 'auto', origin = 'lower'\
    #, vmin = V.min(), vmax = V.max(), cmap = pl.cm.jet)
    
    return cft


def plotBurialHistory_temperature(strat_steps, age_base_all, cellSize_steps, temperature):
    
    #######################################
    # plot burial history
    #######################################
    
    stratMarkers, strat_ages, stratColors = setStratCodeMarkersandColors()
    
    stratMarkers = stratMarkers[:: - 1]
    stratColors = stratColors[:: - 1]
    
    Nsteps = len(strat_steps)
    stratPlotDepths = np.zeros((len(stratMarkers), Nsteps)) ; stratPlotAges = np.zeros((Nsteps))
    
    # go through all tectonic time steps
    for i in xrange(Nsteps):
        stratPlotAges[i] = age_base_all[i]
        
        
        # sample the stratigraphic units present at each time step
        for j, stratName in enumerate(strat_steps[i]):
                                            
            # find the stratigraphic units in the stratigraphic marker list
            for h, stratIndex in enumerate(stratMarkers):
                a = len(stratIndex)
                
                # if stratigrpahic unit is present in marker list: store depth of base of stratigraphic unit
                if stratName[:a] == stratIndex:
                    #print 'plot strat %s' %stratIndex
                    stratPlotDepths[h, i] = cellSize_steps[i, :j].sum()
                elif type(stratIndex) == list:
                    #print 'plot compiled strat:'
                    for stratIndex_item in stratIndex:
                        a = len(stratIndex_item)
                            
                        if stratName[:a] == stratIndex_item:
                            #print stratIndex_item
                            stratPlotDepths[h, i] = cellSize_steps[i][:j].sum()
    
    # plot subsidence:
    for h in xrange(len(stratMarkers) - 1):
        #pl.fill_between(stratPlotAges, (stratPlotDepths[h, :]/1000.), (stratPlotDepths[h + 1, :]/1000.), color = stratColors[h], label = stratMarkers[h])
        leg_burial = pl.plot(stratPlotAges, (stratPlotDepths[h, :]/1000.), color = 'black', lw = 0.5, label = 'Formation depths')
    #pl.plot(stratPlotAges, (stratPlotDepths[h, :]/1000.))
    #pl.fill_between(stratPlotAges, stratPlotDepths[-1, :]/1000., 0, color = stratColors[-1], label = stratMarkers[-1])
    
    print 'plotting temperature grid'
    cft = plotTemperatureGrid(cellSize_steps, age_base_all, temperature)
    
    return cft, leg_burial


def plotProvenanceHistory(AFTsimTime_prov_, AFTsimTemp_prov_, sample_linestyle, sample_colors):

    Nsamples = len(AFTsimTime_prov_)    
    
    geotherm = 30.
    
    for i in xrange(Nsamples):
        Ntimesteps, NprovenanceHistories = np.shape(AFTsimTime_prov_[i])
        
        # draw range of prov. histories
        x0 = AFTsimTime_prov_[i][:, 0].max() - AFTsimTime_prov_[i][:3, 0]
        x1 = AFTsimTime_prov_[i][:,  - 1].max() - AFTsimTime_prov_[i][:3,  - 1]
        
        y_temp0 = AFTsimTemp_prov_[i][:3, 0]
        y_surface0 = y_temp0[2]
        y_burial0 = (y_temp0 - y_surface0)/geotherm
        
        y_temp1 = AFTsimTemp_prov_[i][:3,  - 1]
        y_surface1 = y_temp0[2]
        y_burial1 = (y_temp1 - y_surface1)/geotherm
            
        leg_provFill = pl.fill(np.concatenate((x0, x1[:: - 1])), 
                    np.concatenate((y_burial0, y_burial1[:: - 1])), 
                    color = (0.9, 0.9, 0.9), zorder = 0)
        
        # daw individual end - member scenarios:
        for j in xrange(NprovenanceHistories):
            
            x_ = AFTsimTime_prov_[i][:3, j]
            xmax = AFTsimTime_prov_[i][:, j].max()
            x = xmax - x_
            
            y_temp = AFTsimTemp_prov_[i][:3, j]
            y_surface = y_temp[2]
            y_burial = (y_temp - y_surface)/geotherm
            
            #leg_provPlot = pl.plot(x, y_burial, color = sample_colors[i], ls = sample_linestyle[i], lw = 1.5, zorder = 10)
            leg_provPlot = pl.plot(x, y_burial, color = sample_colors[i], 
                                        ls = '--', lw = 1.5, zorder = 10)
        
    
    return leg_provFill, leg_provPlot


def plotAFTburialHist(cellDepths_present, strat_steps, age_base_all, 
    cellSize_steps,  AFTDepth_base, AFTDepth_top, sample_linestyle, 
    sample_colors):
    
    Nsamples = len(AFTDepth_base) ; Nsteps = len(strat_steps)
    for sampleNo in xrange(Nsamples):
        
        # find sample location index in depth array & stratigraphy list:
        y_start = (AFTDepth_base[sampleNo] + AFTDepth_top[sampleNo])/2.
        AFT_index_top = np.where(y_start > cellDepths_present)[0][-1]
        AFT_index_bottom = np.where(y_start < cellDepths_present)[0][0]
        strat_top = strat_steps[-1][AFT_index_top]
        strat_bottom = strat_steps[-1][AFT_index_bottom]
        
        
        ytop, ybottom = cellDepths_present[AFT_index_top], cellDepths_present[AFT_index_bottom]
        
        factor = (y_start - ytop)/(ybottom - ytop)        
        stratIndex_top, stratIndex_bottom = strat_steps[-1][AFT_index_top], strat_steps[-1][AFT_index_bottom]

        # go through all tectonic time steps and track AFT sample location:
        AFTPlotDepths = np.zeros((2, Nsteps)) ; AFTPlotAges = np.zeros((Nsteps))
        for i in xrange(Nsteps):
            AFTPlotAges[i] = age_base_all[i]
            
            # sample the stratigraphic units present at each time step
            top_index = - 1 ; base_index = - 1
            for j, stratName in enumerate(strat_steps[i]):
                # find the stratigraphic units of the AFT samples
                if stratName[:len(stratIndex_top)] == stratIndex_top:
                    top_index = j
                elif stratName[:len(stratIndex_bottom)] == stratIndex_bottom:
                    base_index = j
                    
            if abs(top_index - base_index) > 1:
                #print 'unconformity in overlying strat unit of AFT sample'
                top_index = base_index
            if top_index!= - 1 and base_index!= - 1:
                AFTPlotDepths[0, i] = cellSize_steps[i, :top_index].sum()
                AFTPlotDepths[1, i] = cellSize_steps[i, :base_index].sum()
            else:
                pass
                #print 'error,  burial depth AFT sample not traceable at step %s of %s' %(i, Nsteps)
                #print 'strat names top: %s,  base: %s' %(stratIndex_top, stratIndex_bottom)
                #exit(1)
        
        # interpolate sample depth from strat unit depths:
        AFTPlotDepths_final = (1 - factor)*AFTPlotDepths[0, :]  +  factor*AFTPlotDepths[1, :]

        # plot AFT sample burial hist:
        pl.plot(AFTPlotAges, AFTPlotDepths_final/1000., 
            color = sample_colors[sampleNo], 
            ls = sample_linestyle[sampleNo], lw = 1.5)
        
    return
    

def plotHeatFlowHistory(age_base_all_, heatflowHist_all_, scenarioLineColors_, scenarioLineStyles_):
    
    pl.plot(age_base_all_, heatflowHist_all_*1000, color = scenarioLineColors_, lw = 2, linestyle = scenarioLineStyles_)

    return


def plotAFTsampleLocations(AFTDepth_base, AFTDepth_top, 
                            sample_linestyle, sample_colors):
    
    ax = pl.gca()
    Nsamples = len(AFTDepth_base)
    xmin, xmax = pl.xlim() ; xdisplay = [] ; ydisplay = []
    xmax = xmax *5
    for sampleNo in xrange(Nsamples):
            
        y = (AFTDepth_base[sampleNo] + AFTDepth_top[sampleNo])/2.
        leg_AFTlocs = pl.fill((xmin, xmin, xmax, xmax )\
        , (AFTDepth_base[sampleNo] ,  AFTDepth_top[sampleNo]\
        ,  AFTDepth_top[sampleNo] ,  AFTDepth_base[sampleNo]),  color='black')
        #leg_AFTlocs = pl.axhline(y, color = sample_colors[sampleNo]\
        #, ls = sample_linestyle[sampleNo], lw = 1.5, zorder = 1)
            
    return leg_AFTlocs

    #leg_AFTlocs=plotAFTsampleLocations(AFTDepth_base, AFTDepth_top, sample_linestyle, sample_colors)

def plotVR(VR_simulated, cellDepths_\
, VRValue, VRValue_std, VRDepth_, residual_VR, lineStyle_, color_, plotGOF=True):
    
    try:
        if VRValue_std.max() > 0:
            leg_vo = pl.errorbar(VRValue, VRDepth_, xerr = VRValue_std, marker = 'v'\
            , ms = 4, linestyle = 'None', color = 'gray', mec = 'black', mfc = 'gray', zorder = 9)
        else:
            leg_vo = pl.errorbar(VRValue, VRDepth_, xerr = 0.1, marker = 'v'\
            , ms = 4, linestyle = 'None', color = 'gray', mec = 'black', mfc = 'gray', zorder = 9)
    except: 
        leg_vo = pl.errorbar((0.2, 1.0), (0.0, 3000.0), visible = False)

    leg_vsim = pl.plot(VR_simulated, cellDepths_, label = 'R0,  sim.', color = color_, linestyle = lineStyle_)

    # show goodness of fit:
    if plotGOF == True:
        tekst = 'GOF=%0.2f' %residual_VR   
        pl.text(0.93, 0.90, tekst, fontsize = 7, bbox=dict(facecolor="white", ec='white', alpha=0.7)\
        , verticalalignment = 'bottom', horizontalalignment = 'right', transform  = pl.gca().transAxes,  zorder=101)
    
    return leg_vo, leg_vsim



def plotTrackLengths(trackLengths_, Nlengths, tracklength_freq,
                    residual, binsize, scenarioLineStyles,
                    scenarioLineColors, verboseNotation=True):

    trackLnHist = pl.hist(trackLengths_, bins = np.arange(0, 18, 1),
                            range = (0, 18), visible = False)
    frequencies_ = trackLnHist[0]
    bins = trackLnHist[1]
  
    # normalize frequencies
    frequencies = (frequencies_.astype(float)/(frequencies_.sum()))
    
    # plot observed track lengths:
    leg_FTlen_obs = pl.bar(bins[:-1], frequencies,
                            facecolor = 'gray', edgecolor = 'darkgray')
    
    #########################################
    # plot simulated AFT length distribution
    #########################################
    a, Nprovenance_ages, N_compositions = np.shape(tracklength_freq)
    MTL_data = trackLengths_.mean()
    MTL_simulated_all = np.zeros((Nprovenance_ages, N_compositions))
    
            
    for i in range(Nprovenance_ages):
        for j in range(N_compositions):
            x = np.arange(0, 20, binsize) - 0.5*binsize
            y = tracklength_freq[:, i, j]*(1./binsize)
            
            #print 'max sim pdf: %s' %y.max()
            #############################################################
            # calculate mean track length from simulated track lenght pdf
            #############################################################
            TL = np.zeros((0))
            for k in range(len(x)):
                for l in range(int(round(y[k]*1000))):
                    TL = np.append(TL, x[k])
            MTL_simulated_all[i, j] = TL.mean()
            
            ########################################
            # plot simulated AFT length distribution
            ########################################            
            pl.plot(x, y, c = scenarioLineColors,
                    ls = scenarioLineStyles, lw = 0.5)

    MTL_simulated_min = MTL_simulated_all[ :, :].min()
    MTL_simulated_max = MTL_simulated_all[ :, :].max()
    
                                
    pl.xlim(0, 18) ; pl.ylim(0, 0.7) ; pl.xticks(np.arange(0, 20, 5))
        
    # determine number of tracks:
    if Nlengths == 0:
        Nlengths = len(trackLengths_)

    # add model fit statistic text    
    tekst = 'N = %.0f'%(Nlengths)
    if verboseNotation == True:
        tekst += ', Mobs = %0.1f'%(trackLengths_.mean())
        tekst += '\nMsim = %0.1f'%(MTL_simulated_min)
        tekst += '- %0.1f'%(MTL_simulated_max)
    ax = pl.gca()
    pl.text( 0.1, 0.94, tekst, fontsize = 7,
                verticalalignment = 'top',
                horizontalalignment = 'left',
                transform  = ax.transAxes )
    
    return leg_FTlen_obs


def plotObservedAFTAges(AFTage_, AFTage_min_, AFTage_max_):
    
    ##########################################################################
    # plot apatite fission track ages. 
    # ages on y axis,  mean 95% confidence intervals on the logarithmic x - axis
    #########################################################################
    
    # sort ages using relative standard error
    SE = ((AFTage_max_ - AFTage_min_)/2.)  # relative standard error
    a = np.argsort(SE)      
    SE = SE[a]
    AFTage = AFTage_[a]
    AFTage_min = AFTage_min_[a]
    AFTage_max = AFTage_max_[a]
    
    # construct error bar plot of ages:
    Nages = len(AFTage)
    ages_range = np.arange(Nages) + 1
    ages_error = np.zeros((2, Nages))
    ages_error[0, :] = - AFTage_min ; ages_error[1, :] = AFTage_max
    
    # plot ages + error bars
    leg_AFTages_obs = pl.errorbar(AFTage, ages_range, xerr = ages_error, ls = 'None', ecolor = 'black', mfc = 'gray', mec = 'black'\
    , marker = 'o', lw = 0.5, ms = 6)
    
    return leg_AFTages_obs


def plotSimulatedAFTages(AFTsample_id_, simulatedAFTages, residual,
                        scenarioLineColors, scenarioLineStyles,
                        ymin=0, ymax=20, verboseNotation=True):
    
    ###################################
    # draw simulated age range:
    ###################################
    Nprovenance_ages, N_compositions = np.shape(simulatedAFTages)
    leg_AFTages_sim_fill = [] ; leg_AFTages_sim_plot = []
        
    # plot box around age range:
    xmin = simulatedAFTages[:, :].min()
    xmax = simulatedAFTages[:, :].max()            
    
    leg_AFTages_sim_fill.append(pl.fill((xmin, xmax, xmax, xmin), 
            (ymin, ymin, ymax, ymax), ec = 'black', fc = 'lightgrey'))
        
    # plot line at individual simulated ages
    for i in xrange(Nprovenance_ages):
        j = 0
        x = simulatedAFTages[i, j]
        x = simulatedAFTages[i, j]
        leg_AFTages_sim_plot.append(pl.plot((x, x), (ymin, ymax),
                            color = scenarioLineColors[0], 
                            ls = scenarioLineStyles[0], lw = 0.5))

    # print text with model fit statistic:
    tekst = '' ; ax = pl.gca()
    
    tekst = 'GOF=%0.2f' %(residual)
    pl.ylim(0, ymax)
    pl.text(0.05, 0.94, tekst, fontsize = 7, va = 'top',  
    transform  = ax.transAxes,  zorder=101)
    
    return leg_AFTages_sim_fill, leg_AFTages_sim_plot


def drawStratigraphicAge(AFTtime_prov_all_, ymax_):

    stratAge = AFTtime_prov_all_[-1, 0] - AFTtime_prov_all_[2, 0]
    leg_stratAge = pl.plot((stratAge, stratAge), (0, ymax_ + 2),
                            color = 'blue', ls = '-', lw = 1.5,
                            zorder=100)

    return leg_stratAge


def removeXtickLabel():
    
    ax = pl.gca()
    pl.setp(ax.get_xticklabels(),  visible = False)

    return


def removeYtickLabel():
    
    ax = pl.gca()
    pl.setp(ax.get_yticklabels(),  visible = False)

    return


def drawUnconformity(y, xmin, xmax, amplitude, freq_, visible_):
    
    Npoints = 50
    
    xdist = xmax - xmin
    x = np.linspace(xmin, xmax, Npoints)
    
    freq = xdist/freq_
    
    
    yrange = np.ones((Npoints))*y
    ywobble = np.ones((Npoints))*y
    
    for i in xrange(Npoints):
        
        dist_ = (math.pi*6/Npoints)*i
        ywobble[i] = (math.sin(dist_))*amplitude
    
    yrange = ywobble + y    
    
    leg_unc = pl.plot(x, yrange, color = 'black', visible = visible_)
    
    #verts_ = []
    #for i in xrange(len(x)):
    #    verts_.append((x[i], yrange[i]))
    #leg_unc = pl.scatter(( - xmax,  - xmin), ( - y,  - y), marker = (verts_, 0), color = 'black')
    #leg_unc = ['~']
    return leg_unc


def getUnconformitySymbol():
    
    xu = np.concatenate((np.linspace( - 1, 1, 20), np.linspace(1,  - 1, 20)))
    yu = np.sin(xu*7)*0.5
    verts = zip(xu, yu)
    
    return verts

#verts = getUnconformitySymbol()

def setStratCodeMarkersandColors():
    
    stratCode_markers = ['D'   ,              'ZE',                'RB',                'RN'               , ['ATAL', 'ATRT', 'ATPO'], ['ATW', 'ATB']               \
    , 'SL'               , 'CK'               , 'NL'            , 'NM'        , 'NU'][:: - 1]
    stratColors =      [(0.584, 0.773, 0.780), (0.980, 0.667, 0.686), (0.627, 0.365, 0.651), (0.878, 0.800, 0.890), (0.000, 0.690, 0.941), (0.45, 0.78, 0.95)\
    , (0.176, 0.710, 0.451), (0.906, 0.929, 0.678), (0.890, 0.745, 0.541), (1.000, 0.808, 0.361), (1.000, 0.867, 0.012)][:: - 1]
    
    #stratColors =      ['grey', 'yellow', 'brown', 'darkred', 'red', 'blue', 'darkgreen', 'lightgreen', 'red', 'orange', 'yellow'][:: - 1]
    
    strat_ages = ['Carboniferous', 'Permian', 'L. Triassic', 'M. - U. Triassic', 'L. Jurassic', 'M. Jurassic', 'U. Jurassic - L. Cretaceous', 'U. Cretaceous'\
    , 'Paleocene', 'Eocene', 'Oligocene - present'][:: - 1]
    
    return stratCode_markers, strat_ages, stratColors
    

def plotStratigraphy(stratCodes, depths, xmin, xmax):
    # add stratigraphy
    
    #print stratCodes
    #exit()
    y_unc2 = - 99999
    stratCode_markers, strat_ages, stratColors = setStratCodeMarkersandColors()
    leg_strat = [] ; labels_strat = [] ; y_unc = - 99999
    for i in range(len(stratCode_markers)):
        startIndex = 0 ; endIndex = 0
        for j in range(len(stratCodes)):
            a = len(stratCode_markers[i])
            if type(stratCode_markers[i]) == list:
                #print 'plot compiled strat:'
                for stratIndex_item in stratCode_markers[i]:
                    a = len(stratIndex_item)
                    if stratCodes[j][:a] == stratIndex_item:
                        endIndex = j
                        if startIndex == 0: startIndex = j
            elif stratCode_markers[i] in stratCodes[j][:a]:
                endIndex = j
                if startIndex == 0: startIndex = j
                
        if startIndex == endIndex and endIndex < 3:
            startIndex = 0
        else:
            startIndex -= 1
        if endIndex!= 0:
            if startIndex == 0: ymin = 0
            else: ymin = depths[startIndex]
            ymax = depths[endIndex]
            leg_strat.append(pl.fill((xmin, xmin, xmax, xmax)\
            , (ymin, ymax, ymax, ymin)\
            , fc = stratColors[i], ec = 'k', lw = 0.5))
            
            if len(labels_strat) == 0:
                labels_strat.append(('Stratigraphy:\n' + strat_ages[i]))
            else:
                labels_strat.append(strat_ages[i])
            if stratCode_markers[i] == 'CK':
                y_unc = ymax
            elif stratCode_markers[i] == 'SL' and y_unc == - 99999:
                y_unc = ymin
            
            elif stratCode_markers[i] == 'AT' and y_unc == - 99999:
                y_unc = ymax
            
            # find pre - permian unconformity:
            if stratCode_markers[i] == 'D' and y_unc2 == - 99999:
                y_unc2 = ymin
            
            print '\t%s  = %0.2f  -  %0.2f ' %(stratCode_markers[i], ymin, ymax)
            #print 'x'*5 + '  ' + (ymin, ymax, ymax, ymin)
        else:
            print '\twarning: %s not found' %stratCode_markers[i]
    # add faults:
    if 'FN' in stratCodes:
        print '\tfound fault in strat column'
        pl.plot((xmin, xmax), (depths[stratCodes.index('FN')], depths[stratCodes.index('FN')]), lw = 2.0, color = 'black')
        #exit()
    else:
        pass
        #print '\tno faults in strat sequence'
        #print stratCodes
        
    # add unconformity
    leg_unc = drawUnconformity(y_unc, xmin, xmax, 0.05, 10, True)
    if y_unc2!= - 99999:
        print '\tdrawing pre - Permian unconformity'
        dummy = drawUnconformity(y_unc2, xmin, xmax, 0.05, 10, True)
        
    #exit()
    return leg_strat, labels_strat, leg_unc, y_unc


def plotDummyStratigraphy():
    xmin = - 5 ; xmax = - 6
    # add stratigraphy
    
    stratCode_markers, strat_ages, stratColors = setStratCodeMarkersandColors()
    leg_strat = [] ; labels_strat = [] ; y_unc = - 99999
    for i in range(len(stratCode_markers)):
        ymin = - i*10 + 5 ; ymax = - i*10 + 10
        
        leg_strat.append(pl.fill((xmin, xmin, xmax, xmax)\
        , (ymin, ymax, ymax, ymin)\
        , fc = stratColors[i], ec = 'k', lw = 0.5, visible = True))
        #if len(labels_strat) == 0:
        #    labels_strat.append(('Stratigraphy:\n' + strat_ages[i]))
        #else:
        labels_strat.append(strat_ages[i])
        if stratCode_markers[i] == 'CK':
            y_unc = ymax
        elif stratCode_markers[i] == 'SL' and y_unc == - 99999:
            y_unc = ymin
        
        elif stratCode_markers[i] == 'AT' and y_unc == - 99999:
            y_unc = ymax
        print '%s  = %s  -  %s ' %(stratCode_markers[i], ymin, ymax)
            #print 'x'*5 + '  ' + (ymin, ymax, ymax, ymin)
    
    leg_unc = drawUnconformity(y_unc, xmin, xmax, 0.05, 10, True)

            
    return leg_strat, labels_strat, leg_unc


def plotMap(plotMode, AFTwellNames, AFTwellLabels):

    # set map extent:
    ymin_map = 51.2 ; ymax_map = 51.9
    xmin_map = 4.8 ; xmax_map = 6.0
    #ymin = 50.9 ; ymax = 52.1
    #xmin = 4.8 ; xmax = 6.2
    # set contour levels:
    cntIntervals = np.arange(0, 2750, 250)

    
    ###################
    # set up basemap
    ###################
    print '\tsetting up basemap'
    bm = Basemap(llcrnrlon = xmin_map, llcrnrlat = ymin_map, urcrnrlon = xmax_map, urcrnrlat = ymax_map, \
    resolution = 'h', lon_0 = (xmin_map + xmax_map)/2, lat_0 = (xmin_map + xmax_map)/2, projection = 'tmerc')

    rasterFiles = '../../GIS/Geokaart_2006/grids_latlong/AT_dik.tif'
    rasterLabel = 'Thickness (m)'

    shapeFileNames = ['../../GIS/Geokaart_2006/faults_latlong/AT_fault.shp']
    shapeFileLabels = ['Faults'] ; shapeFileLS = ['- '] ; shapeFileColors = ['darkgray'] ; shapeFileLW = np.ones((1))*0.5

    print '\treading raster file:'
    print rasterFiles
    cf = plotLib.drawGrid_bm(bm, rasterFiles, rasterLabel, cntIntervals)

    # draw faults
    print '\tdrawing faults'
    plotLib.drawShapeFiles_bm(bm, shapeFileNames, shapeFileLabels, shapeFileColors, shapeFileLS, shapeFileLW)
    
    try:
        hdf5filename = "../../hdf5Data/allRVGdata.h5"
        h5file = tables.openFile(hdf5filename,  "r")
    except:
        hdf5filename = "/home/elco/hdf5Data/allRVGdata.h5"
        h5file = tables.openFile(hdf5filename,  "r")
                    
    #####################
    # draw AFT data locs:
    print '\tget AFT well locations'
    x_AFT, y_AFT = find_AFT_locs(h5file, AFTwellNames)
    x_AFT_map, y_AFT_map = bm(x_AFT, y_AFT)
    # plot well locs:
    pl.scatter(x_AFT_map, y_AFT_map, marker = 'o', label = 'Wells', edgecolor = 'black', facecolor = 'lightgrey', zorder = 80)
            
    # determine well label positions:
    x_AFT_offset, y_AFT_offset = x_AFT.copy(), y_AFT.copy()
    if plotMode!= 'all':
        x_AFT_offset = x_AFT - 0.05 ; y_AFT_offset = y_AFT + 0.05
    else:
        x_AFT_offset = x_AFT ; y_AFT_offset = y_AFT + 0.02
        
    for j, well in enumerate(AFTwellNames):
        if plotMode == 'Geotrack':
            if well == 'BKZ-01' or well == 'SPC-01':
                print 'adjusting offset %s' %well
                x_AFT_offset[j] = x_AFT[j] - 0.25 ; y_AFT_offset[j] = y_AFT[j]
            elif well == 'KDK-01':
                print 'adjusting offset %s' %well
                x_AFT_offset[j] = x_AFT[j] + 0.05 ; y_AFT_offset[j] = y_AFT[j]
            elif well == 'SMG-01':
                print 'adjusting offset %s' %well
                x_AFT_offset[j] = x_AFT[j] - 0.3 ; y_AFT_offset[j] = y_AFT[j] - 0.06
        else:
            if well == 'ALM-01' or well == 'SPC-01' or well == 'SMG-01':
                print 'adjusting offset %s' %well
                x_AFT_offset[j] = x_AFT[j] - 0.08
            elif well == 'WWK-01' or well == 'AND - 06':
                print 'adjusting offset %s' %well
                x_AFT_offset[j] = x_AFT[j] + 0.05

            
    # plot well data
    x_AFT_map, y_AFT_map = bm(x_AFT_offset, y_AFT_offset)
    
    for j, well in enumerate(AFTwellNames):
        print '\t', AFTwellLabels[j], ',  ', well
        print '\tadding well label %s for well %s' %(AFTwellLabels[j], well)
        pl.text(x_AFT_map[j], y_AFT_map[j], AFTwellLabels[j], fontsize = 'xx-small', weight = 'bold', bbox = dict(facecolor = 'lightgrey',  alpha = 0.5), zorder = 100)
            
    print '\tdrawing borders'
    bm.drawcountries()
    bm.drawmeridians(np.arange(int(xmin_map), int(xmax_map) + 0.5, 0.5), labels = [1, 0, 0, 1], fontsize = 'xx-small')
    
    if plotMode != 'Geotrack':
        bm.drawparallels(np.arange(int(ymin_map/0.25)*0.25, int(ymax_map/0.25)*0.25 + 0.25, 0.25), labels = [1, 0, 0, 1], fontsize = 'xx-small')
    else:
        bm.drawparallels(np.arange(int(ymin_map/0.5)*0.5, int(ymax_map/0.5)*0.5 + 0.5,  0.5), labels = [1, 0, 0, 1], fontsize = 'xx-small')
    
    return
    
    
def plotModelResults_profiles_v2(well, strat_steps_allScenarios,
        age_base_all_allScenarios, cellSize_steps_allScenarios, 
        heatflowHist_all_allScenarios, 
        VR_simulated_allScenarios, cellDepths_present_mid_allScenarios,
        VRValue, VRValue_std, VRDepth, residual_VR, 
        AFTsample_id, AFTDepth_base, AFTDepth_top, 
        AFTages, AFTage_min, AFTage_max, residual_AFTages, 
        simulatedAFTages, 
        AFTtime_prov_all, AFTtemp_prov_all,
        trackLengths, trackLengthPDF, residual_AFTlengths, 
        Nplots, plotNo, 
        inversion, heatflow, 
        AFTwellNames, 
        wellStratCode, wellStratDepth_baseTVD, wellStratDepth_topTVD,
        sampleTemp,
        BHTDepth, BHTValue,
        DSTgaugeDepth, DSTintervalDepthTop, DSTintervalDepthBottom,
        DSTtemperature, 
        logDepth, logTemp, BHTcDepth, BHTcTemp, 
        columnNo, y_start_AFT, y_top_strat,
        addLegend=True, addMap=True, maxDepth = 3.5, NlegendCols=1,
        samplePerColumn=5, Ncols = 3, addStratLegend = True):
    
    print '-'*10
    print well
    
    AFTDepth = (AFTDepth_base+AFTDepth_top) /2000.0
    
    scenarioLineStyles = ['solid']
    scenarioLineColors = ['black']
    
    NAFTsamples = len(AFTages)
    
    # set figure margins:
    left = 0.02 ; right = 0.99
    bottom = 0.01 ; top = 1.0 
    # number of columns:
    
    # horizontal sizes of subplots:
    xsize_total = (right-left)/Ncols
    xu = 0.02 * 3/Ncols
    xsize_strat = xu ; xsize_T = xu*2 ; xsize_VR = xu*2
    # space separating AFT and strat/VR/temp plots:
    x_sep = xu
    # space separating AFT plots
    x_sep_AFT = xu/2.0
    # set left indent of columns:
    left_indent = xu*1.5
    
    # calculate size left for AFT plots:
    xsize_left = ( xsize_total - xsize_strat - xsize_T 
                    - xsize_VR - x_sep - left_indent )
    xsize_AFTages = xsize_left * 2.0 / 3.0 - x_sep_AFT
    xsize_AFTlengths = xsize_left * 1.0 / 3.0 - x_sep_AFT
    
    # max number of AFT samples per column
    
    
    # determine vertical size of new subplots:
    ymargin_strat = 0.07 ; ymargin_AFT = 0.01
   
    # separator for lowest y labels:
    y_sep = 0.04

    # vertical size of AFT subplots:
    ysize_AFT = ((top-bottom-y_sep)/samplePerColumn) - ymargin_AFT*2 
    ysize_strat = ysize_AFT * 1.5
    
    # determine total vertical size of well subplots
    ysize_total_AFT = NAFTsamples * (ysize_AFT + ymargin_AFT) + y_sep
    ysize_total_strat = ysize_strat + ymargin_strat + y_sep
    
    # check size & skip to new column if new subplot does not fit:
    y_end_AFT = y_start_AFT - ysize_total_AFT 
    y_end_strat = y_top_strat - ysize_total_strat 
    print 'bottom AFT & strat plots %s, %s' %(y_end_AFT,y_end_strat)
   
    # shift to new column if bottom subplot<bottom figure:
    if y_end_AFT< 0 or y_end_strat < 0:
        y_start_AFT = top
        y_top_strat = top
        y_end_AFT = y_start_AFT - ysize_total_AFT
        y_end_strat = y_top_strat - ysize_total_strat
        columnNo +=1
        
    # set horizontal and vertical position of subplots:
    xplot = xsize_total * columnNo + left_indent
    yplot = y_start_AFT - ysize_AFT    
    # set strat plot at midpoint of AFT samples
    yplot_strat = (y_start_AFT - (NAFTsamples*ysize_AFT) / 2.0 -
                    ysize_strat/2.0)
    
    if yplot_strat+ysize_total_strat > y_top_strat:
        diff = y_top_strat - (yplot_strat+ysize_strat+ymargin_strat)
        if y_top_strat < 1.0:
            diff -= y_sep
        yplot_strat += diff
        
    y_end_strat = yplot_strat - y_sep
    
    print '\tpositioning strat plot at ',yplot_strat
    print 'y_top = %s,  column  = %s' %(y_start_AFT, columnNo)
    print 'start plotting at x = %s,  y = %s for AFT and %s for strat'\
            %(xplot, yplot, yplot_strat)
    
    #####################
    # plot stratigraphy
    #####################
    print 'plot stratigraphy'
    axv = pl.axes((xplot + left, yplot_strat, xsize_strat, ysize_strat))    
    # add fault in well HSW-01
    if well == 'HSW-01':
        for i, s in enumerate(strat_steps_allScenarios[0][-1]):
            if s[:4] == 'ATWD':
                strat_steps_allScenarios[0][-1][i] = 'FN'
                print '\tchanged strat %s to FN' %s  
            else:
                pass
                #print 'no faults in strat %s' %s[:3]
    leg_strat, labels_strat, leg_unc, unc_depth = plotStratigraphy(strat_steps_allScenarios[0][-1]\
    , cellDepths_present_mid_allScenarios[0]/1000.0, 0, 100)

    # draw AFT sample locs:
    leg_AFTlocs=plotAFTsampleLocations(AFTDepth,AFTDepth, ['-']*10, ['darkgray']*10)

    pl.xlim(0, 100)
    pl.ylabel('Depth (km)') 
    #pl.title('(%s) Strat'%(notation[Nsubplot]))
    pl.ylim(maxDepth, 0) ; pl.xticks([])
    
    ###################################
    # plot present day temperature data
    ###################################
    print 'plot temperature data'
    #(xplot, yplot, xsize_strat, ysize_strat)
    axt = pl.axes((xplot+left+xsize_strat, yplot_strat, xsize_T, 
                    ysize_strat))
    leg_BHT_, leg_BHTc_, leg_DST_ = plotTemperatureData(
                                        BHTDepth/1000.0, BHTValue,
                                        DSTgaugeDepth/1000.0,
                                        DSTintervalDepthTop/1000.0,
                                        DSTintervalDepthBottom/1000.0,
                                        DSTtemperature,
                                        logDepth/1000.0, logTemp,
                                        BHTcDepth/1000.0, BHTcTemp)

    # plot simulated temperature data
    Nsamples = len(np.nonzero(sampleTemp[-1, :])[0])
    surfaceTemp = 10.0
    x = np.concatenate((np.array([surfaceTemp]), sampleTemp[-1, :Nsamples]))
    y = np.concatenate((np.array([0]), cellDepths_present_mid_allScenarios[0]/1000.0))
    pl.plot(x, y, color = 'black')    

    # draw AFT sample locs:
    leg_AFTlocs=plotAFTsampleLocations(AFTDepth, AFTDepth, ['-']*10, ['darkgray']*10)
    
    pl.xlim(10.0, 130.0)
    pl.xticks([20, 60, 100], rotation = 45)
    axt.xaxis.grid(True)
    pl.ylim(maxDepth, 0) ; pl.yticks([])
    degree_symbol  = unichr(176) ; labeltext = 'T (%sC)' %degree_symbol
    pl.xlabel(labeltext) #; pl.ylabel('Depth (km)')
    
    # add title for each well
    az = string.ascii_uppercase
    tekst = '('+az[plotNo]+') '
    tekst += well
    tekst += '\nInv=%0.0f,HF=%0.0f' %(inversion, heatflow)
    #tekst += r'$mW m^{ - 2}$'
    pl.title(tekst,  fontsize = 'x-small')
    
    #####################
    # plot VR data
    #####################
    axv = pl.axes((xplot + left + xsize_strat + xsize_T, yplot_strat, xsize_VR, ysize_strat))
    # draw AFT sample locs:
    leg_AFTlocs=plotAFTsampleLocations(AFTDepth, AFTDepth, ['-']*10, ['darkgray']*10)
    if len(VRValue) > 0:
        plotVR_switch = True
        print 'plot VR'
        # draw VR data
        leg_vo, leg_vsim = plotVR( VR_simulated_allScenarios[0],
                        cellDepths_present_mid_allScenarios[0]/1000.,
                        VRValue, VRValue_std, VRDepth/1000.,
                        residual_VR,
                        scenarioLineStyles[0], scenarioLineColors[0] )
        pl.xlabel(r'VR (R$_o$)')
        xvmax=VRValue.max()*1.2
        pl.ylim(maxDepth, 0)
        if xvmax <= 2.0:
            pl.xticks([0.2, 0.5, 1.0, 1.5, 2.0], rotation = 45)
        else:
            pl.xticks([0.2, 1.0, 2.0, 3.0, 4.0])
        axv.xaxis.grid(True)
        if xvmax < 1.1:
            xvmax = 1.1
        pl.xlim(0.0, xvmax)
        removeYtickLabel()
    else:
        plotVR_switch = False
        xa, xb = pl.xlim()
        VRValue = ( - 2000.0,  - 2000.0)
        VRDepth_ = ( - 2000.0,  - 2000.0)
        leg_vo = pl.errorbar(VRValue, VRDepth_, xerr = 0.1, zorder=10)
        scenarioNo = 0
        leg_vsim = pl.plot(VRValue, VRDepth_, label = 'R0,  sim.',
                            color = scenarioLineColors[0],
                            linestyle = scenarioLineStyles[0],
                            zorder=11)
        pl.ylim(maxDepth, 0) ; pl.xlim(xa, xb)
        removeYtickLabel() ; pl.xticks([])
        
    #################
    # plot AFT data:
    #################
    print '\tplot AFT data'
    Nsamples = len(AFTages) ; ax_ages = [] ; con = []
    
    # sort samples from top to bottom:
    AFT_depth_index = np.argsort(AFTDepth)
    
    xplot_AFTages = xplot+left+xsize_strat+xsize_T+xsize_VR+x_sep
    xplot_AFTlengths = xplot+left+xsize_strat+xsize_T+xsize_VR+x_sep+xsize_AFTages
    ymargin = 0.01
    
    for sampleNo in xrange(Nsamples):
        
        print '\t\tsample %i,  %s' %(sampleNo + 1,  AFTsample_id[sampleNo]) 
        
        # set axes
        ax_ages.append(pl.axes((xplot_AFTages, yplot, xsize_AFTages, ysize_AFT - ymargin)))
        
        # set axis limits
        Nages = len(AFTages[sampleNo])
        if Nages < 20: ymax = 20
        else: ymax = Nages
        
        # plot observed AFT ages:
        leg_AFTages_obs = plotObservedAFTAges(AFTages[sampleNo], AFTage_min[sampleNo], AFTage_max[sampleNo])
        
        # plot simulated AFT ages:
        leg_AFTages_sim_fill,  leg_AFTages_sim_plot = plotSimulatedAFTages(AFTsample_id[sampleNo]\
        ,  simulatedAFTages[sampleNo, :, :],  residual_AFTages[sampleNo],  scenarioLineColors, scenarioLineStyles\
        ,  ymax=Nages,  verboseNotation=False)
        
        # draw strat age:
        leg_stratAge = drawStratigraphicAge(AFTtime_prov_all[sampleNo], Nages)
        pl.axvline(x=0, color='black')
        pl.xlim(500, -25)
        if Nages<15:
            pl.ylim(0, 20)
        else:
            pl.ylim(0, Nages+7)
        ax_ages[-1].grid(True)
        pl.yticks([]) ; pl.setp(pl.gca().get_yticklines(),  visible = False)
        
        if sampleNo == (NAFTsamples - 1):
            pl.xlabel('AFT age (Ma)')
        else:
            ax_ages[-1].set_xticklabels([])
        
        ##################################################
        # connect figures to sample locs in VR plot
        ##################################################
        yA = 1 - (AFTDepth[sampleNo])/((maxDepth))
        xyA_ = (1.0, yA) ; xyB_ = (0, 0)
        #xy = (0.5, 1)
        con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_,
                coordsA = "axes fraction", coordsB = "axes_fraction",
                axesA = axv,  axesB = ax_ages[-1], color = 'black'))
        ax_ages[-1].add_artist(con[-1])
        #
        yA = 1 - (AFTDepth[sampleNo])/((maxDepth))
        xyA_ = (1.0, yA) ; xyB_ = (0, 1.0)
        con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_,
                                        coordsA = "axes fraction",
                                        coordsB = "axes_fraction",
                                        axesA=axv,  axesB=ax_ages[-1],
                                        color = 'black'))
        ax_ages[-1].add_artist(con[-1])
        
        ##########################################
        # plot observed and simulated AFT lengths
        ##########################################
        axl = pl.axes((xplot_AFTlengths, yplot, xsize_AFTlengths,
                        ysize_AFT - ymargin))
        #ax = pl.axes((left + width1, y_plot, dx2, dy))
        binsize = 0.25
        # calculate amount of observed length data
        NtrackLn = len(trackLengths[sampleNo])
        if well == 'WWK-01':
            sample = AFTsample_id[sampleNo]
            WWK_samples = ['WWK01-1', 'WWK01-3', 'WWK01-4', 'WWK01-5',
                            'WWK01-6', 'WWK01-7', 'WWK01-8', 'WWK01-9']
            WWK_Ntracks = np.array([8, 0, 101, 87, 12, 22, 47, 12])
            NtrackLn = WWK_Ntracks[WWK_samples.index(sample)]
        else: NtrackLn = 0
            
        # plot trackLengths
        leg_AFTlen_obs = plotTrackLengths(
                            trackLengths[sampleNo], NtrackLn,
                            trackLengthPDF[sampleNo, :, :, :],
                            residual_AFTlengths[sampleNo], binsize,
                            scenarioLineStyles[0],
                            scenarioLineColors[0],
                            verboseNotation=False)
    
        # add axis labels  
        #pl.yticks([0.2, 0.4, 0.6])  ; pl.ylabel('rel. freq.')
        #if sampleNo == 0: pl.title('(%s) Apatite fission\ntrack lengths'%(notation[3])) 
        #if sampleNo == Nsamples - 1:
        #    mu_symbol  = '$\mu' ; labeltext = 'Track length (%sm)' %mu_symbol
        #    pl.xlabel(r'Track Length ($\mu m$)')
        axl.grid(True)
        pl.xticks([5, 10, 15])
        if sampleNo == (NAFTsamples - 1):
            pl.xlabel(r'TL ($\mu m$)')
        else:
            removeXtickLabel()
        pl.yticks([])
        
        ##################################################
        # connect figures to sample locs in AFT ages plot
        ##################################################
        xyA_ = (1.0, 1.0) ; xyB_ = (0, 1.0)
        con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_,  coordsA = "axes fraction"\
        ,  coordsB = "axes_fraction", axesA = ax_ages[-1],  axesB = axl, color = 'black'))
        axl.add_artist(con[-1])
        
        xyA_ = (1.0, 0.0) ; xyB_ = (0, 0.0)
        con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_,  coordsA = "axes fraction"\
        ,  coordsB = "axes_fraction", axesA = ax_ages[-1],  axesB = axl, color = 'black'))
        axl.add_artist(con[-1])
        
        #tekst = 'sample %s' %(sampleNo + 1)
        #pl.text(1.2, 0.2, tekst, fontsize = 'small', rotation = 90, transform  = ax.transAxes)
    
        yplot -= ysize_AFT
        
    
    #####################################
    # add figure with AFT well locations:
    #####################################
    if addMap == True:
        try:
            xplot = xsize_total*columnNo
            yplot = bottom + ((1 - bottom) - ysize_AFT*rowNo)
        
            print 'adding map with well locations'
            tekst = '('+['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'][plotNo+1]+') '
            tekst += 'Well locations'
            AFTwellLabels = AFTwellNames
            
            pl.axes((xplot , yplot-(3.5*ysize_AFT),  xsize_total,  (2.5*ysize_AFT)),  aspect = 'equal')
            plotMap('Geotrack',  AFTwellNames,  AFTwellLabels)
            fontsize_default = 'x-small'
            pl.title(tekst,  fontsize = fontsize_default)
            
            #pdb.set_trace()
        except:
            print 'x' * 10
            print 'warning, failed to add map of sample locations'
        
    if addLegend == True:
        
        ######################################################
        # add dummy strat column with all strat units included
        ######################################################
        x1, x2 = pl.xlim() ; y1, y2 = pl.ylim()
        leg_strat, labels_strat, leg_unc = plotDummyStratigraphy()
        
        leg_fault = pl.plot(( -10,  -9), ( -10,  -9), lw = 2.0, color = 'black')
        leg_header = pl.plot(( -10,  -9), ( -10,  -9), visible = False)
        verts = getUnconformitySymbol()
        #leg_unc = pl.scatter(( - 10,  - 9), ( - 10,  - 9), marker = (verts, 0), s = 30, color = 'black', visible = False)
        leg_unc = leg_header[0]
        leg_header = leg_header[0]

        leg_BHT = pl.scatter(( -10,  -9), ( -10,  -9),  marker = 'o', c = 'lightgrey', s = 3,  label = 'BHT')
        leg_BHTc = pl.scatter(( -10,  -9), ( -10,  -9), marker = 'o', c = 'darkgrey',  label = 'corrected BHT')
        leg_DST = pl.scatter(( -10,  -9), ( -10,  -9), c = 'k',  marker = 's', label = 'DST')

        #leg_BHT = pl.scatter(( 1, 1), ( 1, 1),  marker = 'o', c = 'lightgrey', s = 3,  label = 'BHT')
        #leg_BHTc = pl.scatter(( 1, 1), ( 1, 1), marker = 'o', c = 'darkgrey',  label = 'corrected BHT')
        #leg_DST = pl.scatter(( 1, 1), ( 1, 1), c = 'k',  marker = 's', label = 'DST')


        pl.xlim(x1, x2) ; pl.ylim(y1, y2)
        
        ###############
        # add legend:
        ###############
        print 'creating legend'
        
        # construct legen entries:
        entries = [leg_header]
        if addStratLegend == True:
            for i in xrange(len(leg_strat)):
                entries.append(leg_strat[i])
        entries += [leg_unc, leg_fault, leg_stratAge]
        entries += [leg_header, leg_BHT, leg_BHTc, leg_DST, leg_vo[0], leg_AFTages_obs[0], leg_AFTlen_obs[0]]
        entries += [leg_header, leg_vsim[0], leg_AFTages_sim_plot[0], leg_AFTages_sim_fill[0]]
        
        #except:
        #    print 'error,  failed to create legend entries'
        
        # set labels:
        if addStratLegend == True:
            leg_labels = ['Stratigraphy:'] + labels_strat 
        else:
            leg_labels = []
            
        leg_labels += ['Major unconformity'\
        , 'Fault', 'Strat. Age', 'Observed data:'\
        , 'BHT'\
        , 'Corrected BHT', 'DST temp.'\
        , 'VR reflectance', 'AFT ages', 'AFT lengths'\
        , 'Modeled data:', 'Temperature/VR/Track lengths'\
        , 'AFT age end-members', 'AFT age range']
            
        pl.figlegend(entries, leg_labels, loc = 'lower right', 
                        ncol = NlegendCols)
        
    #pdb.set_trace()
    
    return columnNo, y_end_AFT, y_end_strat
    
    
def plotTemperatureData(BHTDepth, BHTValue, DSTgaugeDepth, DSTintervalDepthTop, DSTintervalDepthBottom, \
    DSTtemperature, logDepth, logTemp, BHTcDepth, BHTcTemp):
    
    ###############################
    #plot temperature data of well
    ###################################
    
    # uncorrected bottom hole temperatures
    try: leg_BHT = pl.errorbar(BHTValue,  BHTDepth,  marker = 'o', ms = 3, c = 'lightgrey', linestyle = 'None', label = 'BHT')
    except: leg_BHT = []
    
    # log temperatures
    try:
        pl.plot(logTemp, logDepth, color = 'darkgray', label = 'log')
    except: pass
        
    # corrected bottom hole temperatures
    try:
        Ndata = len(BHTcTemp) ; x_err = np.ones((Ndata))*5 
        leg_BHTc = pl.errorbar(BHTcTemp,  BHTcDepth, xerr = x_err,  marker = 'o', c = 'darkgrey', linestyle = 'None', label = 'corrected BHT')
    except:
        leg_BHTc = []
    
    # DST's
    try:
        Ndata = len(DSTintervalDepthTop) 
        ybest = (DSTintervalDepthTop + DSTintervalDepthBottom)/2
        y_err = np.zeros((Ndata,  2))
        y_err[:,  0] = ybest - DSTgaugeDepth ; y_err[:,  1] = DSTintervalDepthBottom - ybest
        x = DSTtemperature ; x_err = np.ones((Ndata))*3        
        leg_DST = pl.errorbar(x,  ybest, [y_err[:, 0], y_err[:, 1]], x_err, c = 'k',  linestyle = 'None',  marker = 's', label = 'DST')
        #pl.errorbar(DSTintervalDepthTop, DSTtemperature, xerr = 5, yerr = color = 'black', label = 'dst')
    except:
        leg_DST = []
    
    return leg_BHT, leg_BHTc, leg_DST


def plotModelResults(well, strat_steps_allScenarios, 
    age_base_all_allScenarios, cellSize_steps_allScenarios, 
    heatflowHist_all_allScenarios, cellDepths_present_base_allScenarios, 
    cellDepths_present_mid_allScenarios, temperature, 
    VR_simulated_allScenarios, VRValue, VRValue_std, 
    VRDepth, residual_VR, AFTsample_id, AFTDepth_base, 
    AFTDepth_top, AFTages, AFTage_min, AFTage_max, residual_AFTages, 
    simulatedAFTages, AFTtime_prov_all, AFTtemp_prov_all, 
    trackLengths, trackLengthPDF, residual_AFTlengths, 
    Nplots, plotNo, 
    addHeatFlow = True,  addTrackLengths = True,  
    addProvenance = True,  addSampleNotation = True, 
    addSimulatedAges = True):

    scenarioLineStyles = ['solid']
    scenarioLineColors = ['black']
    
    
    # settings for presentation figure:
    #addHeatFlow = False
    #addTrackLengths = False
    #addProvenance = False
    #addSimulatedAges = False

    Nscenarios = len(age_base_all_allScenarios)
    ysize_all = (1.0/Nplots)
    ybottom_all = 1.0 - ysize_all*(plotNo + 1)
    ybottom_plots = 0.35*ysize_all + ybottom_all
    left = 0.08
    ysize_burial = 0.4/Nplots
    ysize_HF = 0.10/Nplots
    ysize_FT = 1.0/Nplots
    
    if Nplots > 1:
        firstChar = ['A', 'B', 'C', 'D', 'E', 'F'][plotNo]
        notation = []
        for i in xrange(1, 10): notation.append((firstChar + str(i)))
    else:
        notation = ['A', 'B', 'C', 'D', 'E', 'F']
        firstChar = ''
    
    #sample_linestyle = ['-', '-.', ':', '-.', '-']*2
    #sample_colors = ['darkgray', 'black', 'black', 'black', 'black']*2
    sample_linestyle = ['-', '-', '-', '-', '-']*2
    sample_colors = ['black', 'black', 'black', 'black', 'black']*2
    
    #sample_linestyle, sample_colors
    
    print 'creating figure %s/%s of model results' %(plotNo, Nplots)
    print '\tbottom plots: %s ,  %s' %(ybottom_all, ybottom_plots)
    #####################
    # plot burial history
    #####################
    xsize_burial = 0.3
    if len(VRValue)<=1:
        xsize_burial += 0.10
    axb = pl.axes((left, ybottom_plots + ysize_HF + 0.01,
                    xsize_burial, ysize_burial))

    cft, leg_burial = plotBurialHistory_temperature(
            strat_steps_allScenarios[0], age_base_all_allScenarios[0],
            cellSize_steps_allScenarios[0], temperature)
    a, maxDepth = pl.ylim()
    
    # add AFT sample locs:
    pl.scatter(np.zeros((len(AFTDepth_base))), 
            (AFTDepth_base/1000.+AFTDepth_top/1000.)/2.0,  marker='o', 
            color='black')
    
    minTime, maxTime = pl.xlim() ; maxTime = maxTime*1.1
    minTime = -10
    pl.xlim(maxTime, minTime)
        
    #####################################
    # plot AFT sample provenance history 
    #####################################
    if addProvenance == True:
        leg_provFill,  leg_provPlot =\
            plotProvenanceHistory(AFTtime_prov_all, 
                AFTtemp_prov_all, sample_linestyle, sample_colors)
        # find maximum provenance age
        a, maxTime = pl.xlim() ; maxTime = maxTime*1.1
        if maxTime < 450: maxTime = 450
    
    ##################################
    # plot AFT sample burial history
    ##################################
    plotAFTburialHist(cellDepths_present_base_allScenarios[0], 
        strat_steps_allScenarios[0], age_base_all_allScenarios[0], 
        cellSize_steps_allScenarios[0],  AFTDepth_base, AFTDepth_top, 
        sample_linestyle, sample_colors)
    pl.yticks([ -0.5,  0,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5])
    pl.xlim(maxTime, minTime)#; removeXtickLabel()
    pl.ylim(maxDepth,  - 0.25)
    pl.ylabel('Depth (km)') ; pl.title('(%s) Thermal history'%(notation[0])) 
    
    #################
    # plot heat flow
    #################
    if addHeatFlow == True:
        axh = pl.axes((left, ybottom_plots, 0.3, ysize_HF))
        plotHeatFlowHistory(age_base_all_allScenarios[0], 
                            heatflowHist_all_allScenarios[0],
                            scenarioLineColors[0],
                            scenarioLineStyles[0] )
        #
        a, b = pl.ylim()
        if b - a < 30:
            ystep = 5
        elif b - a < 50:
            ystep = 10
        else:
            ystep = 25
        ymin_ = int(a/ystep)*ystep
        ymax_ = int(math.ceil(b/ystep)*ystep)
        pl.yticks(np.arange(ymin_, ymax_ + ystep, ystep))
        pl.xlim(maxTime, 0)
        pl.ylabel(('Heat flow\n' + r'$(mW m^2)$'))
        pl.xlabel('Time (Ma)')
        
    ###############################
    # plot temperature color legend
    ###############################
    if Nplots == 1 or plotNo < Nplots - 1:
        axcb = pl.axes((left, ybottom_plots - 0.08, 0.2, 0.01))
        pl.xticks([]) ; pl.yticks([])
        cb = pl.colorbar(cft, cax = axcb, ax = axb, orientation = 'horizontal')
        degree_symbol  = unichr(176) ; labeltext = 'Temperature (%sC)' %degree_symbol
        cb.set_label(labeltext, fontsize = 'xx-small')
    
    #####################
    # plot VR data
    #####################
    if len(VRValue)>1:
        axv = pl.axes((left + 0.31, ybottom_plots + ysize_HF + 0.01, 0.10, ysize_burial))
        pl.yticks([])
        leg_vo = plotVR(
                        VR_simulated_allScenarios[0],
                        cellDepths_present_mid_allScenarios[0]/1000.,
                        VRValue, VRValue_std, VRDepth/1000.,
                        residual_VR,
                        scenarioLineStyles[0], scenarioLineColors[0])
        pl.title('(%s) VR'%(notation[1])) ; pl.xlabel(r'$R_o$')
        pl.xticks([0, 0.5, 1.0, 1.5, 2.0]) ; xvmin, xvmax = pl.xlim() ; 
        # plot AFT sample locs
        try:
            leg_AFTlocs = plotAFTsampleLocations(AFTDepth_base/1000., 
                    AFTDepth_top/1000., sample_linestyle, sample_colors)
        except:
            print 'failed to plot AFT depths'
            print AFTDepth_base
            print AFTDepth_base/1000.
        
        pl.ylabel('')
        pl.ylim(maxDepth,  - 0.25) ; pl.xlim(xvmin,  xvmax)
    else:
        axv = axb
        
    #################
    # plot AFT data:
    #################
    Nsamples = len(AFTages) ; ax_ages = []
    con = []
    
    # sort samples from top to bottom:
    AFT_depth_index = np.argsort(AFTDepth_base)
    
    for sampleNo in xrange(Nsamples):
        
        Nsamples_ = Nsamples
        if Nsamples_ < 2: Nsamples_ = 2

        
        # define plot locations:
        top = 0.9*ysize_all + ybottom_all
        bottom = 0.08*ysize_all + ybottom_all
        left = 0.55 ; right = 0.97 ; wspace = 0.06 ; hspace = 0.03
        height = (top - bottom)/Nsamples_
        width1 = (right - left)*2./3.
        width2 = (right - left) - width1 
        y_plot = bottom + (Nsamples_ - sampleNo - 1)*height + hspace
        dy = height - hspace
        dx1 = width1 - wspace ; dx2 = width2
        
        # main plot
        ax_ages.append(pl.axes((left + 0.05, y_plot, dx1 - 0.05, dy), 
                                                    frame_on = False))
        
        ############################
        ymax_ = len(AFTages[sampleNo])
        if ymax_ < 20: ymax = 20
                
        ##########################
        # plot observed AFT ages:
        ##########################
        leg_AFTages_obs = plotObservedAFTAges(AFTages[sampleNo], 
                            AFTage_min[sampleNo], AFTage_max[sampleNo])
        
        ##########################
        # plot simulated AFT ages:
        ##########################
        leg_AFTages_sim_fill, leg_AFTages_sim_plot = \
            plotSimulatedAFTages(AFTsample_id[sampleNo], 
                    simulatedAFTages[sampleNo, :, :], 
                    residual_AFTages[sampleNo], scenarioLineColors, 
                    scenarioLineStyles)
        
        # draw strat age:
        leg_stratAge =\
            drawStratigraphicAge(AFTtime_prov_all[sampleNo], ymax_)
        
        # adjust x - axis
        
        pl.xlim(500, 0) ; pl.yticks([])
        pl.setp(pl.gca().get_yticklines(),  visible = False)
        if sampleNo < Nsamples - 1: removeXtickLabel()
        else:
            major_formatter  = pl.FormatStrFormatter('%0.0f')
            pl.gca().xaxis.set_major_formatter(major_formatter)
        pl.ylim(0, ymax_ + 8) ; pl.axhline(0, color = 'black', lw = 1)
        pl.axhline(ymax_ + 8, color = 'black', lw = 1)
        pl.axvline(0, color = 'black', lw = 1) 
        if sampleNo == 0:
            tekst = '(%s) AFT ages'%(notation[2])
            pl.title(tekst)
            
        # main plot
        ax_ages.append(pl.axes((left, y_plot, 0.05, dy), frame_on = False))
        
        ##########################
        # plot observed AFT ages:
        ##########################
        dummy = plotObservedAFTAges(AFTages[sampleNo], 
                            AFTage_min[sampleNo], AFTage_max[sampleNo])
        
        pl.xlim(1000, 500) ; pl.xticks((1000, 500))
        
        if sampleNo < Nsamples - 1: removeXtickLabel()
        else:
            major_formatter  = pl.FormatStrFormatter('%0.0f')
            pl.gca().xaxis.set_major_formatter(major_formatter)

        yticks_ = np.arange(0, (math.ceil(ymax_/10.)*10) + 10, 10)
        pl.yticks(yticks_)  ; pl.ylim(0, ymax_ + 8)
        pl.setp(pl.gca().get_yticklines(),  visible = False)
        pl.axhline(0, color = 'black', lw = 1)
        pl.axhline(ymax_ + 8, color = 'black', lw = 1)
        pl.axvline(1000, color = 'black', lw = 1) 
        
        pl.ylabel('N grains')
        if sampleNo == Nsamples - 1: pl.xlabel('AFT age (Ma)')
        
        ##################################################
        # connect figures to sample locs in VR plot
        ##################################################
        yA = 1 - (AFTDepth_base[sampleNo] + 250)/((maxDepth + 0.25)*1000.0)
        xyA_ = (1.0, yA) ; xyB_ = (0, 0)
        #xy = (0.5, 1)
        con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_, 
                coordsA = "axes fraction",  coordsB = "axes_fraction", 
                axesA = axv,  axesB = ax_ages[-1],  color = 'black', 
                alpha = 0.5))
        ax_ages[-1].add_artist(con[-1])
        
        yA = 1 - (AFTDepth_top[sampleNo] + 250)/((maxDepth + 0.25)*1000.)
        xyA_ = (1.0, yA) ; xyB_ = (0, 1.0)
        con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_, 
                coordsA = "axes fraction",  coordsB = "axes_fraction", 
                axesA = axv,  axesB = ax_ages[-1], color = 'black', 
                alpha = 0.5))
        ax_ages[-1].add_artist(con[-1])
        
        ##########################################
        # plot observed and simulated AFT lengths
        ##########################################
        if addTrackLengths == True:
            ax = pl.axes((left + width1, y_plot, dx2, dy))
            binsize = 0.25
            # calculate amount of observed length data
            NtrackLn = len(trackLengths[sampleNo])
            if well == 'WWK-01':
                sample = AFTsample_id[sampleNo]
                WWK_samples = ['WWK01-1', 'WWK01-3', 'WWK01-4', \
                'WWK01-5', 'WWK01-6', 'WWK01-7', 'WWK01-8', 'WWK01-9']
                WWK_Ntracks = np.array([8, 0, 101, 87, 12, 22, 47, 12])
                NtrackLn = WWK_Ntracks[WWK_samples.index(sample)]
            else: NtrackLn = 0
                
            # plot trackLengths
            leg_AFTlen_obs = plotTrackLengths(trackLengths[sampleNo], 
                        NtrackLn, trackLengthPDF[sampleNo, :, :, :], 
                        residual_AFTlengths[sampleNo], binsize, 
                        scenarioLineStyles[0], scenarioLineColors[0])
        
            # add axis labels  
            pl.yticks([0.2, 0.4, 0.6])  ; pl.ylabel('rel. freq.')
            if sampleNo == 0:
                pl.title('(%s) AFT lengths'%(notation[3])) 
            if sampleNo == Nsamples - 1:
                mu_symbol  = '$\mu'
                labeltext = 'Track length (%sm)' %mu_symbol
                pl.xlabel(r'Track Length ($\mu m$)')
            else: removeXtickLabel()
            
            ##################################################
            # connect figures to sample locs in AFT ages plot
            ##################################################
            xyA_ = (1.0, 1.0) ; xyB_ = (0, 1.0)
            con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_, 
                    coordsA = "axes fraction",  coordsB = "axes_fraction", 
                    axesA = ax_ages[-1],  axesB = ax, color = 'black', 
                    alpha = 0.5))
            ax.add_artist(con[-1])
            
            xyA_ = (1.0, 0.0) ; xyB_ = (0, 0.0)
            con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_, 
                    coordsA = "axes fraction",  coordsB = "axes_fraction", 
                    axesA = ax_ages[-1],  axesB = ax, color = 'black', 
                    alpha = 0.5))
            ax.add_artist(con[-1])
            
            #tekst = 'sample %s' %(sampleNo + 1)
            #pl.text(1.2, 0.2, tekst, fontsize = 'small', rotation = 90, transform  = ax.transAxes)
            
    ###############    
    # add legend
    ###############
    try:
        print 'creating legend'
        global entries
        entries = [leg_burial, leg_provPlot, leg_provFill, leg_AFTlocs, \
                    leg_vo[0][0], leg_AFTages_obs[0], leg_AFTlen_obs[0], \
                    leg_stratAge]
        leg_labels = ['Burial depth', 'Provenance burial scenarios', \
            'Provenance burial depth range', 'Burial depth AFT samples', \
            'VR reflectance', 'Apatite fission track ages', \
            'Apatite fission track lengths', 'Stratigraphic age']
        entries.append(leg_AFTages_sim_plot[0])
        tekst = 'Simulated VR/\ntrack length data'
        leg_labels.append(tekst)
        entries.append(leg_AFTages_sim_fill[0])
        tekst = 'Simulated AFT age range'
        leg_labels.append(tekst)    
    except:
        print 'failed creating legend'
    try:
        pl.figlegend(entries,  leg_labels,  loc = 'lower left',  ncol = 2)
    except:
        print 'warning ,  could not create figure legend'
        
    if Nplots > 1:
        tekst = '(%s) Model results,  well %s'%(firstChar, well)
    else:
        tekst = 'Model results,  well ' + well
    pl.figtext(0.08, ybottom_all + ysize_all - 0.01, tekst, 
                    fontsize = 'x-small', va = 'top')
    
    print 'done creating figure'
    
    return


def plotModelResults_presentation(well, strat_steps_allScenarios, 
    age_base_all_allScenarios, cellSize_steps_allScenarios, 
    heatflowHist_all_allScenarios, cellDepths_present_base_allScenarios, 
    cellDepths_present_mid_allScenarios, temperature, 
    VR_simulated_allScenarios, VRValue, VRValue_std, 
    VRDepth, residual_VR, AFTsample_id, AFTDepth_base, 
    AFTDepth_top, AFTages, AFTage_min, AFTage_max, residual_AFTages, 
    simulatedAFTages, AFTtime_prov_all, AFTtemp_prov_all, 
    trackLengths, trackLengthPDF, residual_AFTlengths, 
    scenarioLineStyles, scenarioLineColors, Nplots, plotNo, 
    addHeatFlow = True,  addTrackLengths = True,  
    addProvenance = True,  addSampleNotation = True, 
    addSimulatedAges = True):

    # settings for presentation figure:
    addHeatFlow = False
    addTrackLengths = False
    addProvenance = False
    addSimulatedAges = False
    
    #Nscenarios = len(age_base_all_allScenarios)
    Nscenarios = 1 ; Nplots = 1 ; plotNo = 0
    ysize_all = (1.0/Nplots)
    ybottom_all = 1.0 - ysize_all*(plotNo + 1)
    ybottom_plots = 0.35*ysize_all + ybottom_all
    left = 0.08
    ysize_burial = 0.4/Nplots
    ysize_HF = 0.10/Nplots
    ysize_FT = 1.0/Nplots
    
    if Nplots > 1:
        firstChar = ['A', 'B', 'C', 'D', 'E', 'F'][plotNo]
        notation = []
        for i in xrange(1, 10): notation.append((firstChar + str(i)))
    else:
        notation = ['A', 'B', 'C', 'D', 'E', 'F']
        firstChar = ''
    
    #sample_linestyle = ['-', '-.', ':', '-.', '-']*2
    #sample_colors = ['darkgray', 'black', 'black', 'black', 'black']*2
    sample_linestyle = ['-', '-', '-', '-', '-']*2
    sample_colors = ['black', 'black', 'black', 'black', 'black']*2
    
    #sample_linestyle, sample_colors
    
    print 'start creating figure %s/%s of model results' %(plotNo, Nplots)
    print '\tbottom plots: %s ,  %s' %(ybottom_all, ybottom_plots)
    #####################
    # plot burial history
    #####################
    xsize_burial = 0.3
    if len(VRValue)<=1:
        xsize_burial += 0.10
    axb = pl.axes((left, ybottom_plots + ysize_HF + 0.01, xsize_burial, ysize_burial))
    cft, leg_burial = plotBurialHistory_temperature(strat_steps_allScenarios[0], age_base_all_allScenarios[0]\
    , cellSize_steps_allScenarios[0], temperature)
    a, maxDepth = pl.ylim()
    
    # add AFT sample locs:
    leg_AFTlocs = pl.scatter(np.zeros((len(AFTDepth_base))), 
            (AFTDepth_base/1000.+AFTDepth_top/1000.)/2.0,  marker='o', 
            color='black')
    
    minTime, maxTime = pl.xlim() ; maxTime = maxTime*1.1
    minTime = -10
    pl.xlim(maxTime, minTime)
        
    #####################################
    # plot AFT sample provenance history 
    #####################################
    if addProvenance == True:
        leg_provFill,  leg_provPlot =\
            plotProvenanceHistory(AFTtime_prov_all, 
                AFTtemp_prov_all, sample_linestyle, sample_colors)
        # find maximum provenance age
        a, maxTime = pl.xlim() ; maxTime = maxTime*1.1
        if maxTime < 450: maxTime = 450
    
    ##################################
    # plot AFT sample burial history
    ##################################
    plotAFTburialHist(cellDepths_present_base_allScenarios[0], 
        strat_steps_allScenarios[0], age_base_all_allScenarios[0], 
        cellSize_steps_allScenarios[0],  AFTDepth_base, AFTDepth_top, 
        sample_linestyle, sample_colors)
    pl.yticks([ -0.5,  0,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5])
    pl.xlim(maxTime, minTime)#; removeXtickLabel()
    pl.ylim(maxDepth,  - 0.25)
    pl.ylabel('Depth (km)')
    pl.title('Burial and thermal history') 
    
        
    ###############################
    # plot temperature color legend
    ###############################
    if Nplots == 1 or plotNo < Nplots - 1:
        
        axcb = pl.axes((left, ybottom_plots,  0.2,  0.04))
        pl.xticks([]) ; pl.yticks([])
        cb = pl.colorbar(cft, cax = axcb, ax = axb,  
                        orientation = 'horizontal')
        degree_symbol  = unichr(176)
        labeltext = 'Temperature (%sC)' %degree_symbol
        cb.set_label(labeltext, fontsize = 'xx-small')
    
    #####################
    # plot VR data
    #####################
    if len(VRValue)>1:
        axv = pl.axes((left + 0.31, ybottom_plots + ysize_HF + 0.01, 0.10, ysize_burial))
        pl.yticks([])
        leg_vo = plotVR(VR_simulated_allScenarios[0], cellDepths_present_mid_allScenarios[0]/1000.\
        , VRValue, VRValue_std, VRDepth/1000., residual_VR[0]\
        , scenarioLineStyles[0], scenarioLineColors[0])
        pl.title('VR\nReflectance')
        pl.xlabel(r'$R_o$')
        pl.xticks([0, 0.5, 1.0, 1.5, 2.0])
        xvmin, xvmax = pl.xlim() ; 
        # plot AFT sample locs
        plotAFTsampleLocations(AFTDepth_base/1000., 
                    AFTDepth_top/1000., sample_linestyle, sample_colors)
        
        pl.ylabel('')
        pl.ylim(maxDepth,  - 0.25) ; pl.xlim(xvmin,  xvmax)
    else:
        axv = axb
        
    #################
    # plot AFT data:
    #################
    Nsamples = len(AFTages) ; ax_ages = []
    con = []
    
    # sort samples from top to bottom:
    AFT_depth_index = np.argsort(AFTDepth_base)
    
    for sampleNo in xrange(Nsamples):
        
        Nsamples_ = Nsamples
        if Nsamples_ < 2: Nsamples_ = 2

        
        # define plot locations:
        top = 0.9*ysize_all + ybottom_all
        bottom = 0.08*ysize_all + ybottom_all
        left = 0.55 ; right = 0.97 ; wspace = 0.06 ; hspace = 0.03
        height = (top - bottom)/Nsamples_
        width1 = (right - left)*2./3.
        width2 = (right - left) - width1 
        y_plot = bottom + (Nsamples_ - sampleNo - 1)*height + hspace
        dy = height - hspace
        dx1 = width1 - wspace ; dx2 = width2
        
        # main plot
        ax_ages.append(pl.axes((left + 0.05, y_plot, dx1 - 0.05, dy), 
                                                    frame_on = False))
        
        ############################
        ymax_ = len(AFTages[sampleNo])
        if ymax_ < 20: ymax = 20
                
        ##########################
        # plot observed AFT ages:
        ##########################
        leg_AFTages_obs = plotObservedAFTAges(AFTages[sampleNo], 
                            AFTage_min[sampleNo], AFTage_max[sampleNo])
        
        ##########################
        # plot simulated AFT ages:
        ##########################
        leg_AFTages_sim_fill, leg_AFTages_sim_plot = \
            plotSimulatedAFTages(AFTsample_id[sampleNo], 
                    simulatedAFTages[:, sampleNo, :, :], 
                    residual_AFTages[:, sampleNo], scenarioLineColors, 
                    scenarioLineStyles)
        
        # draw strat age:
        leg_stratAge =\
            drawStratigraphicAge(AFTtime_prov_all[sampleNo], ymax_)
        
        # adjust x - axis
    
        # set gridlines
        ax_ages[-1].xaxis.grid(True)
        
        pl.xlim(500, 0) ; pl.yticks([])
        pl.setp(pl.gca().get_yticklines(),  visible = False)
        if sampleNo < Nsamples - 1: removeXtickLabel()
        else:
            major_formatter  = pl.FormatStrFormatter('%0.0f')
            pl.gca().xaxis.set_major_formatter(major_formatter)
        pl.ylim(0, ymax_ + 8) ; pl.axhline(0, color = 'black', lw = 1)
        pl.axhline(ymax_ + 8, color = 'black', lw = 1)
        pl.axvline(0, color = 'black', lw = 1) 
        if sampleNo == 0:
            tekst = 'Apatite fission\ntrack ages'
            pl.title(tekst)
        
        # main plot
        ax_ages.append(pl.axes((left, y_plot, 0.05, dy), frame_on = False))
        
        ##########################
        # plot observed AFT ages:
        ##########################
        dummy = plotObservedAFTAges(AFTages[sampleNo], 
                            AFTage_min[sampleNo], AFTage_max[sampleNo])
        
        pl.xlim(1000, 500) ; pl.xticks((1000, 500))
        
        if sampleNo < Nsamples - 1: removeXtickLabel()
        else:
            major_formatter  = pl.FormatStrFormatter('%0.0f')
            pl.gca().xaxis.set_major_formatter(major_formatter)

        yticks_ = np.arange(0, (math.ceil(ymax_/10.)*10) + 10, 10)
        pl.yticks(yticks_)  ; pl.ylim(0, ymax_ + 8)
        pl.setp(pl.gca().get_yticklines(),  visible = False)
        pl.axhline(0, color = 'black', lw = 1)
        pl.axhline(ymax_ + 8, color = 'black', lw = 1)
        pl.axvline(1000, color = 'black', lw = 1) 
        
        # set gridlines
        ax_ages[-1].xaxis.grid(True)
        
        pl.ylabel('N grains')
        if sampleNo == Nsamples - 1: pl.xlabel('AFT age (Ma)')
        
        ##################################################
        # connect figures to sample locs in VR plot
        ##################################################
        yA = 1 - (AFTDepth_base[sampleNo] + 250)/((maxDepth + 0.25)*1000.0)
        xyA_ = (1.0, yA) ; xyB_ = (0, 0)
        #xy = (0.5, 1)
        con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_, 
                coordsA = "axes fraction",  coordsB = "axes_fraction", 
                axesA = axv,  axesB = ax_ages[-1],  color = 'black', 
                alpha = 0.5))
        ax_ages[-1].add_artist(con[-1])
        
        yA = 1 - (AFTDepth_top[sampleNo] + 250)/((maxDepth + 0.25)*1000.)
        xyA_ = (1.0, yA) ; xyB_ = (0, 1.0)
        con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_, 
                coordsA = "axes fraction",  coordsB = "axes_fraction", 
                axesA = axv,  axesB = ax_ages[-1], color = 'black', 
                alpha = 0.5))
        ax_ages[-1].add_artist(con[-1])
        
        ##########################################
        # plot observed and simulated AFT lengths
        ##########################################
        if addTrackLengths == True:
            ax = pl.axes((left + width1, y_plot, dx2, dy))
            binsize = 0.25
            # calculate amount of observed length data
            NtrackLn = len(trackLengths[sampleNo])
            if well == 'WWK-01':
                sample = AFTsample_id[sampleNo]
                WWK_samples = ['WWK01 - 1', 'WWK01 - 3', 'WWK01 - 4', \
                'WWK01 - 5', 'WWK01 - 6', 'WWK01 - 7', 'WWK01 - 8', 'WWK01 - 9']
                WWK_Ntracks = np.array([8, 0, 101, 87, 12, 22, 47, 12])
                NtrackLn = WWK_Ntracks[WWK_samples.index(sample)]
            else: NtrackLn = 0
                
            # plot trackLengths
            leg_AFTlen_obs = plotTrackLengths(trackLengths[sampleNo], 
                        NtrackLn, trackLengthPDF[:, sampleNo, :, :, :], 
                        residual_AFTlengths[:, sampleNo], binsize, 
                        scenarioLineStyles[0], scenarioLineColors[0])
            
            # add axis labels  
            pl.yticks([0.2, 0.4, 0.6])  ; pl.ylabel('rel. freq.')
            if sampleNo == 0:
                pl.title('(%s) Apatite fission\ntrack lengths'%(notation[3])) 
            if sampleNo == Nsamples - 1:
                mu_symbol  = '$\mu'
                labeltext = 'Track length (%sm)' %mu_symbol
                pl.xlabel(r'Track Length ($\mu m$)')
            else: removeXtickLabel()
            
            ##################################################
            # connect figures to sample locs in AFT ages plot
            ##################################################
            xyA_ = (1.0, 1.0) ; xyB_ = (0, 1.0)
            con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_, 
                    coordsA = "axes fraction",  coordsB = "axes_fraction", 
                    axesA = ax_ages[-1],  axesB = ax, color = 'black', 
                    alpha = 0.5))
            ax.add_artist(con[-1])
            
            xyA_ = (1.0, 0.0) ; xyB_ = (0, 0.0)
            con.append(patches.ConnectionPatch(xyA = xyA_,  xyB = xyB_, 
                    coordsA = "axes fraction",  coordsB = "axes_fraction", 
                    axesA = ax_ages[-1],  axesB = ax, color = 'black', 
                    alpha = 0.5))
            ax.add_artist(con[-1])
            
            #tekst = 'sample %s' %(sampleNo + 1)
            #pl.text(1.2, 0.2, tekst, fontsize = 'small', rotation = 90, transform  = ax.transAxes)
            
        #
        
    ###############    
    # add legend
    ###############
    print 'creating legend'
    global entries
    entries = [leg_burial,  leg_AFTlocs, \
                leg_AFTages_obs[0], \
                leg_stratAge]
    leg_labels = ['Burial depth', 'Depth AFT samples', \
        'Apatite fission track ages', \
        'Stratigraphic age']
    entries.append(leg_AFTages_sim_plot[0])
    tekst = 'Simulated VR/\nAFT data'
    leg_labels.append(tekst)
    entries.append(leg_AFTages_sim_fill[0])
    tekst = 'Simulated AFT age range'
    leg_labels.append(tekst)    
    
    try:
        pl.figlegend(entries,  leg_labels,  loc = 'lower left',  ncol = 2)
    except:
        print 'warning ,  could not create figure legend'
        
    if Nplots > 1:
        tekst = '(%s) Model results,  well %s'%(firstChar, well)
    else:
        tekst = 'Model results,  well ' + well
    pl.figtext(0.08, ybottom_all + ysize_all - 0.01, tekst, 
                    fontsize = 'x-small', va = 'top')
    
    print 'done creating figure'
    
    return




def parameterCrossPlot(wellNames, parameter1, parameter2, GOF_all,
                        GOF_second, GOF_T, wells_plot, maxExhumation,
                        title_, 
                        arangeSubplots=True,  createColorBar=False,
                        xmin=0.1,  xmax=0.5,  ymin=0.08,
                        convertToKm=True):
        
    # find unique well entries
    #wells_plot = np.unique(np.asarray(wellNames))

    Nwells = len(wells_plot)
    if Nwells < 9:
        Nrows = math.ceil(Nwells/2.0) ; Ncols = 2
    else:
        Ncols = math.floor(math.sqrt(Nwells))
        Nrows = math.ceil(Nwells/Ncols)
        
    if Nrows < 1: Nrows = 1 ; Ncols = 1
    print 'rows:', Nrows, ' cols: ', Ncols

    # specify contour levels:
    V = np.arange( - 0.1, 1.1, 0.1)

    leg_entries = [] ; labels = []
        
    for i, well in enumerate(wells_plot):
        if arangeSubplots == True:
            axp = pl.subplot(Nrows, Ncols, i + 1)
        
        condition_a = np.asarray(wellNames) == well

        # select well data:
        wellIndex = np.where(condition_a == True)[0]
        
        # find grid cell size
        dummy = np.unique(parameter1)
        step1 = dummy[1] - dummy[0]
        
        if len(wellIndex) > 0:
            
            # construct grid values for contouring:
            xi  = np.arange(parameter1[wellIndex].min(),
                            parameter1[wellIndex].max() + step1, step1)
            dummy = np.unique(parameter2)
            step2 = dummy[1] - dummy[0]
            yi  = np.arange(parameter2[wellIndex].min(),
                            parameter2[wellIndex].max() + step2, step2)
            zi  = matplotlib.mlab.griddata( parameter1[wellIndex],
                                            parameter2[wellIndex],
                                            GOF_all[wellIndex], xi, yi)
            # plot contours:
            cf  = pl.contourf(xi, yi, zi, V, cmap = pl.cm.jet_r)
        
            # plot boundary acceptable and non-acceptable fits
            cf2 = pl.contour(   xi, yi, zi, np.array([ - 0.7, 0.7]),
                                color = 'black', lw = 1.5)
            
            # plot line of max. exhumation estimate:
            leg_ex = pl.axvline( x=maxExhumation[i],lw=2,color='black')
            
            # plot lines at confines of GOF temperature:
            if GOF_T[wellIndex].max() > 0.7:
                T_index = wellIndex[np.where(GOF_T[wellIndex]>0.7)[0]]
                HFmax = parameter2[T_index].max()
                HFmin = parameter2[T_index].min()
                leg_HFmax = pl.axhline( y=HFmax, lw=2, ls=':',
                                        color='black')
                leg_HFmin = pl.axhline( y=HFmin, lw=2, ls='--',
                                        color='black')
                
            #plot location of best fit
            best_fit_index = wellIndex[np.argmax(GOF_all[wellIndex])]
            leg_bestFit = pl.scatter( parameter1[best_fit_index],
                        parameter2[best_fit_index],
                        color='blue', edgecolor='black',s=12 )
            
        # plot data points
        leg_data = pl.scatter(  parameter1, parameter2,
                                color = 'gray', s = 1)
        
        # convert x-axis label to km
        if convertToKm == True:
            locs, labels = pl.xticks()
            labels = (locs/1000.0).tolist()
            pl.xticks(locs,labels)
            pl.xlabel('Exhumation (km)')
        else:
            pl.xlabel('Exhumation (m)')
            
        pl.ylabel(r'Heat flow (mW m$^{ - 2}$)')
        a, b = pl.xlim() ; pl.xlim(0, b)
        pl.xlim(parameter1.min(), parameter1.max())
        pl.ylim(parameter2.min(), parameter2.max())
        
        header_ = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
        title2_ = '(%s) well %s' %(header_[i], well)
        pl.title(title2_)
        
        # create colorbar
        if arangeSubplots == True:
            if i >= (len(wells_plot) - 1):
                print 'drawing colorbar'
                axl_test = pl.subplot(Nrows, Ncols, i + 2, frameon = False)
                pl.xticks([]) ; pl.yticks([])
                ax_dimension = axl_test.get_position().get_points()
                xmin, ymin = ax_dimension[0]
                xmax, ymax = ax_dimension[1]
                ydist = 0.015
                xdist = xmax - xmin
                xdist = xdist/6*5
                axl = pl.axes((xmin, (ymin + ymax)/2.0, xdist, ydist))
                cb = pl.colorbar(cf, ax = axp, cax = axl, orientation = 'horizontal', format = '%0.1f')
                #cb.set_label('Mean GOF,  apatite fission track age\nand VR reflectance data.')
                
                label = 'GOF'
                cb.set_label(label)
        
        elif createColorBar == True:
            print 'drawing colorbar'
            ydist = 0.02
            xdist = xmax - xmin
            axp = pl.gca()
            axl = pl.axes((xmin,  ymin , xdist,  ydist))
            cb = pl.colorbar(cf, ax = axp, cax = axl, orientation = 'horizontal', format = '%0.1f')
            #cb.set_label('Mean GOF,  apatite fission track age\nand VR reflectance data.')
            cb.set_label('GOF')
    
        # create legend
        if well == wells_plot[-1]:
            print 'create legend'
            labels = [  'parameter set', 'best fit',
                        'maximum exhumation',
                        'min. present-day heat flow',
                        'max. present-day heat flow']
            entries = [leg_data,leg_bestFit,leg_ex,leg_HFmin, leg_HFmax]
            pl.figlegend(entries,labels,loc='lower right')
    
    # 
    xl = pl.xlim() ; yl = pl.ylim()
    #dummy = pl.plot(( - 500,  - 600), ( - 500,  - 600), color = 'black')
    pl.xlim(xl) ; pl.ylim(yl)
    #leg_entries.append(dummy)
    #tekst = 'Numbered contours = \n' + r'mean error apatite fission track lengths ($\mu m$)'
    #labels.append(r'Numbered contours = \nmean error apatite fission track lengths ($\mu m$)')
    #labels.append(tekst)
    if arangeSubplots == True:
        pl.subplots_adjust(hspace = 0.45, wspace = 0.35, bottom = 0.05,
                            top = 0.95, left = 0.07, right = 0.95)
    #pl.figlegend(leg_entries, labels, loc = 'lower right', ncol = 1)
    #pl.suptitle(title_)
    
    return
