import sys, math, random, tables, os, pickle, time, pdb

import numpy as np
import pylab as pl
import scipy.optimize as opt
import scipy.stats as stats

try:
    import fipy
except:
    print 'warning, failed to import Fipy library'

# local libraries:
sys.path.append("../Libraries/")

# import apatite fission track libraries:
try:
    import AFTlib, AFTannealingLib
except:
    print 'warning, failed to import AFT libraries'

# import vitrinite reflectacne algorithm
try:
    import easyro
except:
    print 'warning, failed to import easryro module'
global verbose 
verbose = False


def linearInterpolation(targetDepth, sampleDepths_, sampleValues):
    
    # linearInterpolation(999, arange(1000), (arange(1000)+100))
    if len(sampleDepths_) != len(sampleValues):
        print 'error'
        print 'unequal length of input arrays linear interpolation'
        print len(sampleDepths_)
        print len(sampleValues)
        sys.exit()
    
    if targetDepth in sampleDepths_:
        targetValue = sampleValues[np.where(targetDepth==sampleDepths_)]

    elif targetDepth >= sampleDepths_[-1]:
        targetValue = sampleValues[-1]
    
    elif targetDepth<sampleDepths_[-1]:
        depthIndex_low = np.where(sampleDepths_<targetDepth)[0][-1]
        depthIndex_high = np.where(sampleDepths_>targetDepth)[0][0]

        fraction = ((sampleDepths_[depthIndex_high]-targetDepth)/
        (sampleDepths_[depthIndex_high]-sampleDepths_[depthIndex_low]))
        
        if fraction<0 or fraction>1:
            print 'error linear interpolation subroutine'
            print 'fraction is %s' %fraction
            print targetDepth
            print sampleDepths_[depthIndex_high]
            print sampleDepths_[depthIndex_low]
            exit()
        
        targetValue = ( fraction*sampleValues[depthIndex_high] + 
                        (1-fraction)*sampleValues[depthIndex_low] )
    
    try: targetValue = targetValue[0]
    except: pass
    
    return targetValue


def findStrat(wellStratTable, well):
    
    # extract well stratigraphy from hdf5 database    
    
    faultStratNo = [] # reset fault list
    
    condition = "well  ==  '%s'" %well
    wellStratDepth_topTVD = np.asarray( [row['top_TVD']
                        for row in wellStratTable.where(condition)] )
    wellStratDepth_baseTVD = np.asarray( [row['bottom_TVD']
                        for row in wellStratTable.where(condition)] )
    wellStratCode = ( [row['unit'] for row 
                                in wellStratTable.where(condition)] )

    # resample Zechstein units,  no distinction between members:
    ZE_count = 0
    ZE_list = np.array([])
    for i in xrange(len(wellStratCode)):
        if wellStratCode[i][:2] == 'ZE':
            if ZE_count == 0:
                ZE_count = i
                ZE_top =  wellStratDepth_topTVD[i]
                wellStratCode[i] = 'ZE'
            else:
                ZE_list = np.append(ZE_list, i)
            ZE_base = wellStratDepth_baseTVD[i]
    
    if ZE_count != 0:
        print 'condense %s Zechstein members to one formation entry:'\
            %(len(ZE_list)+1)
        wellStratDepth_baseTVD[ZE_count] = ZE_base
        # delete member entries if zechstein members are present:
        for i in xrange(len(wellStratCode)-1, -1, -1):
            if i in ZE_list:
                del wellStratCode[i]
                wellStratDepth_topTVD =\
                        np.delete(wellStratDepth_topTVD, i)
                wellStratDepth_baseTVD =\
                        np.delete(wellStratDepth_baseTVD, i)       

    return wellStratCode, wellStratDepth_baseTVD, wellStratDepth_topTVD


def findStratAges(wellStratCode, stratCodesTable):
    
    # find stratigraphic ages of each stratigraphic interval in a well 
    # (wellStratTable) from hdf5 database of strat. units 

    Nstrat = len(wellStratCode)
    age_base_array = np.zeros((Nstrat), dtype = float)
    age_top_array = np.zeros((Nstrat), dtype = float)
    
    for counter in range(Nstrat):
        Nlen = len(str(wellStratCode[counter]))
        if str(wellStratCode[counter])[0] == '-':
            wellStratCode[counter] = wellStratCode[counter][1:]
        # treat formation subdivisions same as parent formation
        if Nlen>2 and wellStratCode[counter][-2] == 'f':   
            wellStratCode[counter] = wellStratCode[counter][:-3]
        if Nlen>3 and wellStratCode[counter][-4:-2] == '-f':  
            wellStratCode[counter] = wellStratCode[counter][:-4]
        if Nlen>2 and wellStratCode[counter][-2] == '-':   
            wellStratCode[counter] = wellStratCode[counter][:-2]
        # use previous member id for hiatus
        if wellStratCode[counter] == 'hiatus':
            wellStratCode[counter] = wellStratCode[counter-1] 

        code_member = wellStratCode[counter]
        if code_member[:5] == 'FAULT' or code_member[:3] == 'FLT':
            code_member = wellStratCode[counter-1]
        if code_member[:2] == 'ZE': code_member = 'ZE'
        if code_member[:4] == 'SLDN': code_member = 'SLDN'
        
        # hdf5 database search condition
        condition = "code_fm  ==  '%s'" %code_member
        try:
            age_base_array[counter] = [row['age_bottom'] 
                    for row in stratCodesTable.where(condition)][0]
        except:
            print 'warning,  cannot find "%s" in strat table %s'\
                %(condition, stratCodesTable)
            pdb.pm()
        age_top_array[counter] = [row['age_top'] 
                    for row in stratCodesTable.where(condition)][0]
        
    return age_base_array, age_top_array


def findTemperatureData(well, h5file):
    
    """
    extract temperaature data from hdf5 database
    """
    
    # find temperature data tables:
    BHTTable = h5file.root.temperature.raw_bht
    DSTTable = h5file.root.temperature.dst
    logTable = h5file.root.temperature.temp_logs
    BHTcTable = h5file.root.temperature.bht_model_results
    
    condition = "well  ==  '%s'" %well
    
    #find uncorrected BHT data:
    BHTDepth = np.asarray([row['depth_TVD'] 
                                for row in BHTTable.where(condition)])
    BHTValue = np.asarray([row['temperature'] 
                                for row in BHTTable.where(condition)])

    #find DST data:
    DSTgaugeDepth = np.asarray([row['depth_gauge_TVD'] 
                        for row in DSTTable.where(condition)])
    DSTintervalDepthTop = np.asarray([row['top_interval_TVD']
                        for row in DSTTable.where(condition)])
    DSTintervalDepthBottom = np.asarray([row['bottom_interval_TVD']
                        for row in DSTTable.where(condition)])
    DSTtemperature = np.asarray([row['temperature']
                        for row in DSTTable.where(condition)])

    #find temperature log data:    
    logDepth = np.asarray([row['depth_TVD']
                        for row in logTable.where(condition)])
    logTemp = np.asarray([row['temperature']
                        for row in logTable.where(condition)])

    #find corrected BHTs:
    BHTcDepth = np.asarray([row['depth_TVD']
                                for row in BHTcTable.where(condition)])
    BHTcTemp = np.asarray([row['temperature']
                                for row in BHTcTable.where(condition)])
    
    return BHTDepth, BHTValue, DSTgaugeDepth, DSTintervalDepthTop,\
            DSTintervalDepthBottom, DSTtemperature,\
            logDepth, logTemp, BHTcDepth, BHTcTemp



def calculateSteadyStateTemperatures(q_,Tsurface_,depthArray,Kbulk_):
        
    tempArray_ = np.zeros((len(depthArray)))
    for j in range(len(depthArray)):
        if j == 0:
            distance = depthArray[j]
            Kbulk_harm = Kbulk_[j]
            tempArray_[j] = Tsurface_+(q_*distance/Kbulk_harm)
        
        else:
            distance = depthArray[j]-depthArray[j-1]
            Kbulk_harm = stats.hmean(Kbulk_[(j-1):(j+1)])
            tempArray_[j] = tempArray_[j-1]+(q_*distance/Kbulk_harm)
        
    return tempArray_


def add_nonDeposition(subsidence_code, strat_units, 
    cellSize_initial_all, cellSize_present_all, age_top, age_base):
    
    print 'adding non deposition time steps'
    
    j = 1
    while j<(len(subsidence_code) - 1):
        
        if (age_base[j] - age_top[j-1])< -1e-4:
            
            #print 'time lag between %s :%0.2f and %s: %0.2f'\
            # %(strat_units[j-1], age_top[j-1], strat_units[j], age_base[j])
            
            subsidence_code.insert(j, 'non_deposition')
            #print 'add non - deposition phase after strat unit %s - %s' %(j - 1, strat_units[j-1])
            
            # update stratigraphic codes table
            strat_units.insert(j, strat_units[j-1])
            
            # add eroded formations to initial cell size array:
            cellSize_initial_all = np.insert(cellSize_initial_all, j, 0)
            
            #
            cellSize_present_all = np.insert(cellSize_present_all, j, 0)
            
            # add base and top ages:
            newtop = age_base[j]
            newbase = age_top[j-1]
            if newtop >=  newbase:
                print 'error,  <=  zero length non_deposition between %s and %s' %(strat_units[j], strat_units[j-1])
                print 'corrected by setting start age of %s to %s' %(strat_units[j], age_top[j-1])
                age_base[j] = age_top[j-1]            
                #exit()
            age_top = np.insert(age_top, j, newtop)
            age_base = np.insert(age_base, j, newbase)
            j += 1
        j += 1
    #exit()
    #print strat_units
    #exit()
    return subsidence_code, strat_units, cellSize_initial_all, cellSize_present_all, age_top, age_base




def calculateVR(sampleTemperatures, subsidence_code, Nsamples, 
                        age_base, VR_timestep):
    
    print 'sample temperature histories for VR reflectance calculations'
    
    Nsamples = len(np.nonzero(sampleTemperatures[-1, :])[0])
    Ntimesteps, a = np.shape(sampleTemperatures)
    VR_simulated = np.zeros((Nsamples))
    
    colorArray = np.linspace(0.1, 0.8, Nsamples)
    
    simulation_time_all = []
    simulated_temp_all = []
    
    
    for i in xrange(Nsamples):
        
        # sample all temperature values of stratigraphic unit:
        simulated_temp = sampleTemperatures[np.nonzero(sampleTemperatures[:, i]), i]
        simulated_temp = np.insert(simulated_temp, 0, simulated_temp[0])
        
        # find corresponding times:
        simulation_time_ = age_base[np.nonzero(sampleTemperatures[:, i])]
        
        # recalculate time to My after deposition instead of My before present
        maxTime = simulation_time_.max()
        simulation_time_ =  - simulation_time_ + simulation_time_.max()
        simulation_time_ = np.append(simulation_time_, maxTime) # add last time step

        # check for errors:
        duration_test = simulation_time_[1:] - simulation_time_[: - 1]
        if duration_test.min() <=  0.:
            print '!!! error sampling temperature histories in VR function'
            print 'unit no %s,  <=  0 My duration timestep, time step array:' %i
            print simulation_time_
            print 'age array taken at indices:'
            print np.nonzero(sampleTemperatures[:, i])
            exit()
        
        # check time step size
        time_diff = abs(simulation_time_[1:] - simulation_time_[: - 1])
        
        # resample if steps > 0.5 My
        if time_diff.max()>VR_timestep:
            
            # resample in smaller timesteps
            factor = int(np.ceil(time_diff.max()/VR_timestep))
            
            simulation_time_fine = np.zeros(((len(simulation_time_) - 1)*factor + 1))
            simulated_temp_fine = np.zeros(((len(simulated_temp) - 1)*factor + 1))
            simulation_time_fine[-1] = simulation_time_[-1]
            simulated_temp_fine[-1] = simulated_temp[-1]
            
            for j in xrange(len(simulation_time_fine) - 1):
                
                step = j - (int(j/factor)*factor)
                
                firstStep = int(j/factor)
                nextStep = int((j + factor)/factor)
                
                #print 'first, next: %s, %s' %(firstStep, nextStep)
                
                diff_time = (simulation_time_[nextStep] - simulation_time_[firstStep])
                diff_temp = (simulated_temp[nextStep] - simulated_temp[firstStep])
            
                # interpolate new time step
                simulation_time_fine[j] = simulation_time_[firstStep] + step*diff_time/factor
                
                # interpolate new temperature:
                simulated_temp_fine[j] = simulated_temp[firstStep] + step*diff_temp/factor
                
            # calculate VR reflectance
            
            if simulation_time_[-1] !=  simulation_time_fine[-1] or simulated_temp[-1] !=  simulated_temp_fine[-1]:
                print 'warning'
                print simulation_time_
                print simulation_time_fine
            
            simulation_time_all.append(simulation_time_fine)
            simulated_temp_all.append(simulated_temp_fine)
    
        else:
            
            simulation_time_all.append(simulation_time_)
            simulated_temp_all.append(simulated_temp)
    
    print 'calculate VR reflectance:'
    #startTime = time.time()
    #VR_simulated  = basinModLib.easyRoFunc_list_fast(simulation_time_all, simulated_temp_all)
    
    VR_simulated = np.zeros((len(simulation_time_all)))
    for i in xrange(len(simulation_time_all)):
        Nsteps_ = len(simulation_time_all[i])
        if Nsteps_ > 4:
            ro = np.zeros((Nsteps_))
            ro = easyro.easyro(ro, simulation_time_all[i]
                                        , simulated_temp_all[i]
                                        , Nsteps_ )
            VR_simulated[i] = ro[-1].copy()
        else:
            VR_simulated[i] = 0.2
    #print '- total duration: %i sec' %(time.time() - startTime)

    #pdb.set_trace()
        
    return VR_simulated
    




def setGridCells(wellStratDepth_baseTVD):
    # align grid cells with formation boundaries in well
    
    Ncells =  len(wellStratDepth_baseTVD) 
    cellSize_present = np.zeros((Ncells))
    cellSize_present[0] = wellStratDepth_baseTVD[0]
    for cell in xrange(1, Ncells):
        cellSize_present[cell] = wellStratDepth_baseTVD[cell] - wellStratDepth_baseTVD[cell - 1]

    # eliminate 0 cell sizes
    cellSize_present = cellSize_present + np.where(cellSize_present == 0, 1, 0)
    
    return cellSize_present


def addErodedAltenaFms(well, erodedATTable, memberTable, 
                        subsidence_code, strat_units, 
                        cellSize_initial_all, cellSize_present_all,
                        age_base, age_top, lateCret_erosion_start,
                        lateCret_erosion_end, maxErosion):
    
    print 'add eroded Jurassic Altena fm. members'
    
    # add the eroded (Jurassic) Altena group formations/members,  +  the associated later erosional time steps
    # assuming a complete succession was deposited in the well prior to erosion
    # normal succession  +  average member thickness taken from hdf5 database
    
    condition = "well  ==  '%s'" %well
    top_preserved_AT = [row['top_preserved_AT'] for row in erodedATTable.where(condition)][0]
    normalSuccession = memberTable.cols.stratigraphic_code[:].tolist()
    ATcellSize = memberTable.cols.thickness_avg[:]

    ATage_base = np.array([199.6, 183., 177.6, 170.6, 170.1, 168.2, 165.4, 164.7, 164., 162.8, 161.2]) # base ages of AT formation members
    ATage_top = np.array([183., 177.6, 170.6, 170.1, 168.2, 165.4, 164.7, 164., 162.8, 161.2, 159.15]) # top ages of AT formation members

    ATaddIndex = strat_units.index(top_preserved_AT) + 1

    # set number of eroded formations
    N_formations_to_add = len(normalSuccession) - normalSuccession.index(top_preserved_AT) - 1
    print '\tadding %s strat units.' %N_formations_to_add

    # check if top unit much lower thickness than average
    diff_topATunit = ATcellSize[ - N_formations_to_add - 1] - cellSize_initial_all[ATaddIndex - 1]
    if diff_topATunit>ATcellSize[ - N_formations_to_add - 1]*0.25:
        N_formations_to_add += 1
        print '\tadd to thickness of topmost preserved unit:'
        print '\told size: %0.0f' %ATcellSize[ - N_formations_to_add]
        ATcellSize[ - N_formations_to_add] = diff_topATunit
        print '\tnew cell size: %0.0f' %ATcellSize[ - N_formations_to_add]
        
    # set array with start times of erosion of AT  + SLDNA formations
    lateCret_erosion_base = np.arange(lateCret_erosion_start, lateCret_erosion_end, (lateCret_erosion_end - lateCret_erosion_start)/(N_formations_to_add + 1))
    lateCret_erosion_top = np.append(lateCret_erosion_base[1:], lateCret_erosion_end)

    # add eroded AT formation members:
    totalATerosion = 0
    print '\tmax eroded thickness: %0.0f' %(maxErosion)
    for j in xrange(N_formations_to_add):
        if totalATerosion >=  (maxErosion - 1):
            print '\tmax erosion reached. not adding %s' %normalSuccession[ - (j + 1)]

        else:
            print '\tadd unit %s & erosion phase, eroded thickness = %0.0f' %(normalSuccession[ - (j + 1)], totalATerosion)
            
            # update subsidence codes:
            subsidence_code.insert(ATaddIndex, 'add_formation')
            subsidence_code.insert(ATaddIndex + (j*2) + 1, 'erosion')
            
            # update stratigraphic codes table
            strat_units.insert(ATaddIndex, normalSuccession[ - (j + 1)])
            strat_units.insert(ATaddIndex + (j*2) + 1, ('-' + normalSuccession[ - (j + 1)]))
            
            # add eroded formations to initial cell size array:
            cellSize_initial_all = np.insert(cellSize_initial_all, ATaddIndex, ATcellSize[ - (j + 1)])
            cellSize_initial_all = np.insert(cellSize_initial_all, ATaddIndex + (j*2) + 1, - ATcellSize[ - (j + 1)])
            totalATerosion += ATcellSize[ - (j + 1)]
            
            cellSize_present_all = np.insert(cellSize_present_all, ATaddIndex, ATcellSize[ - (j + 1)])
            cellSize_present_all = np.insert(cellSize_present_all, ATaddIndex + (j*2) + 1, - ATcellSize[ - (j + 1)])
                
            # add base ages of Altena formations:
            age_base = np.insert(age_base, ATaddIndex, ATage_base[ - (j + 1)])
            age_base = np.insert(age_base, ATaddIndex + j + 1, lateCret_erosion_base[ - (j + 1)])
            
            # add top ages of Altena formations:
            age_top = np.insert(age_top, ATaddIndex, ATage_top[ - (j + 1)])
            age_top = np.insert(age_top, ATaddIndex + j + 1, lateCret_erosion_top[ - (j + 1)])
        
    return subsidence_code, strat_units, cellSize_initial_all, cellSize_present_all, age_base, age_top\
    , lateCret_erosion_base, lateCret_erosion_top, totalATerosion


    
def addErodedSLDNAfm(subsidence_code, strat_units,
                    cellSize_initial_all, cellSize_present_all, 
                    age_base, age_top, lateCret_erosion_start,
                    lateCret_erosion_end, lateCret_erosion_base,
                    lateCret_erosion_top, lateCret_erosion,
                    totalATerosion):
    
    print 'adding eroded lower Cretaceous Schieland fm.'
    
    # add lower cretaceous SLDNA formation if not present, or add eroded section to reach a user - specified thickness (lateCret_erosion)
    # erosional timespan defined by lateCret_erosion_start, erosion_age_end
    # unless Jurassic Altena (AT) formations included in erosion event, then
    # lateCret_erosion_base, lateCret_erosion_top passed on from AT erosion subroutine
    
    totalSLerosion = 0
    SLDNAsize_new = (lateCret_erosion - totalATerosion) # new cell size of SLDNA fm.
    
    # decompact erosion estimate:
    #phi_SLDNA = porosityDepthFunc_feb2010(SLDNAsize_new, 0, array([0.0, 0.4, 0.3, 0.3, 0, 0, 0, 0, 0, 0, 0, 0]), 'x', False\
    #, Phi0, Phic)
    
    if 'SLDNA' not in strat_units and SLDNAsize_new>0:
        #print 'add SLDNA formation'
        ATaddIndex = strat_units.index('ATBRU') + 1
        j = 0
        totalSLerosion = lateCret_erosion - totalATerosion
        
        # update subsidence codes:
        subsidence_code.insert(ATaddIndex, 'add_formation')
        subsidence_code.insert(ATaddIndex + (j*2) + 1, 'erosion')
        
        # update stratigraphic codes table
        strat_units.insert(ATaddIndex, 'SLDNA')
        strat_units.insert(ATaddIndex + (j*2) + 1, '-SLDNA')
        
        # add eroded formations to initial cell size array:
        cellSize_initial_all = np.insert(cellSize_initial_all, ATaddIndex, SLDNAsize_new)
        cellSize_initial_all = np.insert(cellSize_initial_all, ATaddIndex + (j*2) + 1, - SLDNAsize_new)
        
        #
        cellSize_present_all = np.insert(cellSize_present_all, ATaddIndex, SLDNAsize_new)
        cellSize_present_all = np.insert(cellSize_present_all, ATaddIndex + (j*2) + 1, - SLDNAsize_new)
        
        # add base ages of Altena formations:
        age_base = np.insert(age_base, ATaddIndex, 156.4)
        age_base = np.insert(age_base, ATaddIndex + j + 1, lateCret_erosion_base[0])
        
        # add top ages of Altena formations:
        age_top = np.insert(age_top, ATaddIndex, 97.)
        age_top = np.insert(age_top, ATaddIndex + j + 1, lateCret_erosion_top[0])
    
    elif SLDNAsize_new>0:
    
        #print 'modify SLDNA'
        SLDNAIndex = strat_units.index('SLDNA')
        ATaddIndex = strat_units.index('SLDNA') + 1
        
        SLDNAsize_present_decomp = cellSize_initial_all[SLDNAIndex] # current decompacted thickness SLDNA fm.
        SLDNAsize_present_decomp = cellSize_present_all[SLDNAIndex]
        SLDNAsize_add = SLDNAsize_new
        
        totalSLerosion = 0
        #if SLDNA2size>0:
        if SLDNAsize_add>0:
            
            totalSLerosion = SLDNAsize_add
            SLDNA_erosion_fraction = SLDNAsize_add/(SLDNAsize_present_decomp + SLDNAsize_add)
            
            #print 'present size SLDNA: ', SLDNAsize_present_decomp
            #print 'added section of SLDNA fm.: ', SLDNAsize_add
            #print 'fraction: ', SLDNA_erosion_fraction
            #print 'new age:'
            #print 97. + (157. - 97.)*SLDNA_erosion_fraction
            
            if SLDNA_erosion_fraction<0 or SLDNA_erosion_fraction>1:
                
                sys.exit('error calculating eroded thickness of SLDNA fm.')
               
            j = 0
            
            # update subsidence codes:
            subsidence_code.insert(ATaddIndex, 'add_formation')
            subsidence_code.insert(ATaddIndex + (j*2) + 1, 'erosion')
            
            # update stratigraphic codes table
            strat_units.insert(ATaddIndex, 'SLDNA - 2')
            strat_units.insert(ATaddIndex + (j*2) + 1, '-SLDNA - 2')
            
            # add eroded formations to initial cell size array:
            cellSize_initial_all = np.insert(cellSize_initial_all, ATaddIndex, SLDNAsize_add)
            cellSize_initial_all = np.insert(cellSize_initial_all, ATaddIndex + (j*2) + 1, - SLDNAsize_add)
            
            # add eroded formations to present day cell size array:
            cellSize_present_all = np.insert(cellSize_present_all, ATaddIndex, SLDNAsize_add)
            cellSize_present_all = np.insert(cellSize_present_all, ATaddIndex + (j*2) + 1, - SLDNAsize_add)
            
            # add base and top ages of newly added SLDNA fraction formations:
            age_base = np.insert(age_base, ATaddIndex, 97. + (157. - 97.)*SLDNA_erosion_fraction)
            age_top[ATaddIndex - 1] = (97. + (157. - 97.)*SLDNA_erosion_fraction)
            
            # add top ages of Altena formations:
            age_top = np.insert(age_top, ATaddIndex, 97.)
            
            # add erosion phase
            age_base = np.insert(age_base, ATaddIndex + j + 1, lateCret_erosion_start)
            age_top = np.insert(age_top, ATaddIndex + j + 1, lateCret_erosion_end)

     # check duration:
    duration = age_base - age_top
    if duration.min() <=  0:
        print 'error preparing well subsidence data, <0 duration timestep'
        print SLDNAsize_new
        print SLDNAsize_add
        print SLDNAsize_present_decomp 
        print ATaddIndex
        print strat_units[ATaddIndex - 5:ATaddIndex + 5]
        print age_base[ATaddIndex - 5:ATaddIndex + 5]
        print age_top[ATaddIndex - 5:ATaddIndex + 5]
        print SLDNA_erosion_fraction
        
        exit()

    return subsidence_code, strat_units, cellSize_initial_all, cellSize_present_all, age_base, age_top, totalSLerosion

    

def correctAges(age_base, age_top):
    
    for i in xrange(1, len(age_base)):
        duration_test = age_top[i - 1] - age_base[i]
        if duration_test !=  0.:
            #print 'c'        
            if (age_base[i - 1] - age_base[i])>0:
                age_top[i - 1] = age_base[i]
            else:
                age_top[i - 1] = (age_top[i - 1] + age_base[i])/2.
                age_base[i] = age_top[i - 1]
            
            if verbose:
                print 'age corrected, array pos %s, difference %s' %(i, duration_test)
    
    age_top[-1] = 0.
    
    return age_base, age_top
    


def calculateSimulationTime(age_base_sample, age_top_sample):            
    #simulation_time =  - (age_top_sample.copy()) + age_base_sample[0]
    simulation_time = (age_base_sample.copy() + age_top_sample.copy())/2.0
    #print simulation_time
    #exit()
    return simulation_time



def calculateSimulatedAFT_precise(AFTDepth_base_, AFTDepth_top_,
    sampleTemperatures, age_base, age_top, cellDepths_present,
    surfaceTemp, strat_steps, subsidence_code):
    
    #print 'calculating temperature history fission track samples'
    
    ## find simulated temperature history of apatite fission track samples
    #showSampleTemp(sampleTemperatures)
    
    # find stratigraphy corresponding to aft sample
    d = ((AFTDepth_base_ + AFTDepth_top_)/2)    # AFT sample depth
    depthIndex2 = np.argmin(np.where(cellDepths_present<d, 1, 0))  # corresponding depth in well  +  simulated temperature file 
    
    if d>cellDepths_present[-1]:
        depthIndex2 = len(cellDepths_present) - 1
    
    depthIndex1 = depthIndex2 - 1  # corresponding depth in well  +  simulated temperature file 
    
    if depthIndex1<0:
        print 'error trying to find depth %s' %d
        print 'depthIndex2 = %s' %depthIndex2
        depthIndex1 = 0
        
    depthFrac = (d - cellDepths_present[depthIndex1])/(cellDepths_present[depthIndex2] - cellDepths_present[depthIndex1])
    if depthFrac<0 or depthFrac>1:
        print 'error interpolation of data from different depths when sampling AFT temperature history'
        print 'depth fraction = ', depthFrac
        print 'target depth =  ', d
        print 'cell depth top =  %s, index =  %s' %(cellDepths_present[depthIndex1], depthIndex1)
        print 'cell depth bottom =  %s, index =  %s' %(cellDepths_present[depthIndex2], depthIndex2)
        depthfrac = 0.0
        print 'setting depth to ', depthfrac
       
    
    a = len(strat_steps[-1])
    if depthIndex2 >=  a:
        depthIndex2 = depthIndex1
    
    
    if verbose == True:
        try:
            print 'AFT sample in between strat units %s and %s, depth fractions %0.2f and %0.2f'\
            %(strat_steps[-1][depthIndex1], strat_steps[-1][depthIndex2], (1 - depthFrac), depthFrac)
        except:
            print 'error, cannot find AFT sample'
            print depthIndex2
            print len(strat_steps[-1])
            pdb.set_trace()
            
    # find first timestep in temperature history array:
    Nstart2 = 0 ; Nstart1 = 0
    Nstart1 = np.where(sampleTemperatures[:, depthIndex1]>0)[0][0]
    Nstart2 = np.where(sampleTemperatures[:, depthIndex2]>0)[0][0]
    
    # check if unconformity between two strata, if yes, follow temperature history of lowest unit
    if abs(Nstart1 - Nstart2)>1:
        if verbose:
            print 'unconformity between strat units, following T history lowest unit'
            print 'changing start step from %s to %s' %(Nstart1, Nstart2)
        Nstart1 = Nstart2
        depthIndex1 = depthIndex2
        
    # get temperature history of AFT sample from T history overlying and underlying grid cell
    sampleTemperature_AFT = depthFrac*sampleTemperatures[Nstart2:, depthIndex2].copy()\
     + (1 - depthFrac)*sampleTemperatures[Nstart2:, depthIndex1].copy()
    
    if Nstart1 !=  Nstart2:
        diff = Nstart1 - Nstart2
        
        try:
            sampleTemperature_AFT[:diff] =\
                    sampleTemperatures[Nstart2:Nstart1, depthIndex2]
            sampleTemperature_AFT[diff:] = ( depthFrac * 
                            sampleTemperatures[Nstart1:, depthIndex2] +
                            (1 - depthFrac) *
                            sampleTemperatures[Nstart1:, depthIndex1] )
        except:
            print 'error finding target depth for AFT sample at depth %s' %d
            exit()
    
    if verbose:
        #print '-'*10
        print 'simulated temperatures:'
        print sampleTemperature_AFT
    #exit()
    
    simulation_time = calculateSimulationTime(age_base[Nstart2:], age_top[Nstart2:])
    # add first step at surface temperature:
    #simulation_time = np.insert(simulation_time, 0, 0)
    #sampleTemperature_AFT = np.insert(sampleTemperature_AFT, 0, surfaceTemp)
    
    # check for 0 temperatures:
    if sampleTemperature_AFT[0] == 0:
        while sampleTemperature_AFT[0] == 0:
            sampleTemperature_AFT = sampleTemperature_AFT[1:].copy()
            simulation_time = simulation_time[1:].copy()
    
    if len(simulation_time) !=  len(sampleTemperature_AFT):
        print 'error'
        print len(simulation_time)
        print len(sampleTemperatures_)
        print len(age_top)
        exit()
    
    #d = simulation_time[1:] - simulation_time[: - 1]
    #if d.min() <=  0:
    #    print 'error,  <=  0 timestep duration at step %s' %np.argmin(d)
    #    exit()
        
    return simulation_time, sampleTemperature_AFT



def calculateGOFages(AFTages_, AFTages_min_, AFTages_max_,
                        simulatedAFTages_):
    
    #print 'calculating fit statistic fission track ages'
    
    minAge = simulatedAFTages_[:, :].min() ;maxAge = simulatedAFTages_[:, :].max()
    P_AFTage_ = np.zeros((len(AFTages_)))
            
    for ageDataNo in xrange(len(AFTages_)):
            
            # if any of simulated ages within bounds: set residual to 0: 
            if AFTages_[ageDataNo] >=  minAge and AFTages_[ageDataNo] <=  maxAge:
                difference_norm = 0
            
            # calculate residual
            else:
                if (minAge - AFTages_[ageDataNo]) < (AFTages_[ageDataNo] - maxAge):
                   targetAge = minAge
                else:
                    targetAge = maxAge
                difference = targetAge - AFTages_[ageDataNo]
            
                # normalize difference using 95% perc confidence interval of observed fission track ages
                if difference<0:
                    ageRange = AFTages_min_[ageDataNo]
                else:
                    ageRange = AFTages_max_[ageDataNo]
                difference_norm = abs(difference/ageRange)   
                
            # calculate probability of drawing simulated age from normal distribution
            # z - test probability:
            P_AFTage_[ageDataNo] = (1 - stats.norm.cdf(abs(difference_norm)))*2
    
    return P_AFTage_.mean()
    
    

def calculateAFTagesPDF(bins, AFTages, AFTages_min, AFTages_max):
    
    print 'calculating pdf of fission track ages'
    
    AFTages_min_ = AFTages + AFTages_min
    AFTages_max_ = AFTages + AFTages_max
    
    # calculate gamma value (i.e. normally distributed transform of AFT age equation)
    gamma = AFTlib.calculate_gamma_from_AFTage(AFTages)
    #print 'gamma :', gamma
    
    # correct zero ages:
    zero_index = np.where(AFTages == 0)[0]
    gamma[zero_index] = 1e-10
    
    # calculate standard error of gamma:
    gamma_max = AFTlib.calculate_gamma_from_AFTage(AFTages_max_)
    gamma_std = (gamma_max - gamma)/2.0
    
    #print 'AFT max ages:', AFTages_max_
    #print 'gamma of  + 2SE age ', gamma_max
    #print 'std of gamma: ', gamma_std
    
    # calculate pdf of gamma:
    gamma_pdf = np.zeros((len(bins)))
    for i in xrange(len(AFTages)):
        gamma_pdf += stats.norm.pdf(bins, gamma[i], gamma_std[i])
    
    # normalize pdf to a value of 1:
    gamma_pdf = gamma_pdf/(gamma_pdf.sum())
    
    return gamma_pdf

#calculateGOFages(AFTages_, AFTages_min_, AFTages_max_, simulatedAFTages_)
def calculateResidual_AFTages(AFTages, AFTages_min, AFTages_max,
                                simulatedAFTages, doProb):
    
    print 'calculating fit statistic AFT ages'
    
    ##########################################################
    # calculate residual of simulated vs measured AFT ages:
    ##########################################################

    a, Nprovenance_ages, Ncl = np.shape(simulatedAFTages)  
    NAFTsamples = len(AFTages)
    P_AFTage = np.zeros((NAFTsamples)) ; diff_max = 0
    
    ############################
    # calculate pdf of AFT ages:
    ############################
    
    # set gamma distribution bins:
    binRange = 30 ; binSize = 0.01 ; Nbins = binRange/binSize
    bins = np.arange(0, binRange, binSize) ; bins[0] = 1.0e-10
    gamma_pdf_all = np.zeros((NAFTsamples, Nbins))
    
    # calculate AFT ages of bins:
    age_bins = AFTlib.calculate_AFTage_from_gamma(bins)
    
    residual = np.zeros((NAFTsamples)) ; GOF_prob = np.zeros((NAFTsamples))
    
    colors = ['red', 'blue', 'green', 'brown', 'gray']
    
    for sampleNo in xrange(NAFTsamples):
        
        GOF_prob[sampleNo] = calculateGOFages(AFTages[sampleNo], AFTages_min[sampleNo], AFTages_max[sampleNo], simulatedAFTages[sampleNo])
        
        
        # calculate residual:
        residual[sampleNo] = (1.0 - GOF_prob[sampleNo])*100.
        print '\tresidual = %0.2f' %residual[sampleNo]
        age_difference = 0
        try:
            if simulatedAFTages[sampleNo].min()<AFTages[sampleNo].min():
                age_difference += AFTages[sampleNo].min() - simulatedAFTages[sampleNo].min()
        except:
            print 'error evaluating ages: '
            print 'sample number: ', sampleNo
            print 'AFT ages: ', AFTages
            print 'sim. ages: ', simulatedAFTages
            exit()
        if simulatedAFTages[sampleNo].max()>AFTages[sampleNo].max():
            age_difference += simulatedAFTages[sampleNo].max() - AFTages[sampleNo].max()
        
        residual[sampleNo] += abs(age_difference)
        
        print '\tGOF =  %0.2f' %(GOF_prob[sampleNo])
        if verbose:
            print '\trange of sim ages in pdf: %0.0f - %0.0f' %(simulatedAFTages[sampleNo].min(), simulatedAFTages[sampleNo].max())
            print '\trange of pdf ages: %0.0f - %0.0f' %(AFTages[sampleNo].min(), AFTages[sampleNo].max()) 
            print '\t - >added %0.2f to residual' %age_difference
            print '\tresidual: %0.2f' %residual[sampleNo]
            
    if doProb == True:
        return GOF_prob
    else:
        return residual.mean()



def calculateResidual_AFTlengths(trackLengthPDF, trackLengths, 
                                    doProb, binsize):

    print 'calculating residual AFT length data'
    
    #####################################################################################
    # use Kolmogorov - Smirnov test to determine fit of simulated to measured track lengths:
    #####################################################################################

    # determine number of chloride content scenarios:
    a, bins, Nprovenance_ages, Ncl = np.shape(trackLengthPDF)
    
    NrandomSamples = 1000 # number of samples to create a random number array from probability density function
    NAFTsamples = len(trackLengthPDF)
    P_ks = np.zeros((NAFTsamples, Nprovenance_ages, Ncl))
    
    
    for i in xrange(NAFTsamples):
        for j in xrange(Nprovenance_ages):
            for k in xrange(Ncl):
                
                trackLengthPDF_s = trackLengthPDF[i, :, j, k]
                trackLengths_data = trackLengths[i]
            
                if trackLengthPDF_s.max() !=  0:
                    
                    # generate array with track lengths from simulated probability density function:
                    trackLength_sim = np.zeros((0))
                    for l in xrange(int(20/binsize)):
                        for freq in xrange(int(trackLengthPDF_s[l]*NrandomSamples)):
                            #print 'add %s track lengths in bin %s' %(int(trackLengthPDF_s[l]*NrandomSamples), l)
                            simLn = l*binsize + random.random()
                            trackLength_sim = np.append(trackLength_sim, simLn)
                    
                    
                    # scipy kolmogorov smirnof test on 2 randomly distributed variables:
                    D, P_ks[i, j, k] = stats.ks_2samp(trackLengths_data, trackLength_sim)
            
                else:                
                    print 'warning, no result KS test, reset to 0' ; P_ks[i, j, k] = 0
            
                #print 'AFT lengths: kolmogorov - smirnof probability sample %s = %s' %(i, P_ks[i]) 
                if np.isnan(P_ks[i, j, k]) == True:
                    #print 'warning, no result KS test, reset to 0' ; P_ks[i, j, k] = 0
                    pass
        saveLengthsPlot = False
        if saveLengthsPlot == True:
            pl.clf()
            pl.subplot(NAFTsamples, 2, (i*2) + 1)
            pl.hist(trackLengths_data)
            pl.subplot(NAFTsamples, 2, (i*2) + 2)
            try:
                pl.hist(trackLength_sim)
            except:
                pass
            pl.savefig('trackLn_fit_test.eps')
        
        
        
    # select max probabilities for each simulation:
    P_ks_1D = np.zeros((NAFTsamples))
    for i in xrange(NAFTsamples):
        P_ks_1D[i] = P_ks[i, :, :].max()
        
    # calculate residual ( = value to minimize during calibration)
    
    residual_AFTlengths = 1 - P_ks_1D[:].mean()
    GOF_AFTlengths = P_ks_1D[:]
    
    if doProb == False: return residual_AFTlengths
    else: return GOF_AFTlengths


def calculateDifferenceMeanTrackLengths(trackLengths_observed,
                                        tracklength_freq, binsize):
    
    print 'calculating difference mean track lengths'
    
    ##############################################################################################
    # calculate difference between observed and simulated mean track length and standard deviation
    ##############################################################################################

    
    Nsamples, a, Nprovenance_ages, N_compositions = np.shape(tracklength_freq)
    MTL_diff = np.zeros((Nsamples))
    MTL = np.zeros((Nsamples))
    MTL_observed = np.zeros((Nsamples))
    TLSTD = np.zeros((Nsamples))
    MTL_simulated = np.zeros((Nsamples, Nprovenance_ages, N_compositions))
    MTL_simulated_min = np.zeros((Nsamples))
    MTL_simulated_max = np.zeros((Nsamples))
    
    # go through all apatite fission track samples available for well:
    for h in xrange(Nsamples):
        
        trackLengths_observed_ = trackLengths_observed[h]   # get observed lengths
        MTL_data = trackLengths_observed_.mean()   # calculate mean:
        MTL_observed[h] = MTL_data
        
        # go through all simulated provenance age and apatite composition scenarios:
        MTL_ = 1000
        
        for i in xrange(Nprovenance_ages):
            for j in xrange(N_compositions):
                
                # resample simulated track lengths:
                x = np.arange(0, 20, binsize) - 0.5*binsize
                y = tracklength_freq[h, :, i, j]*(1./binsize)
                
                ###################################################
                # calculate simulated mean track length from track lenght pdf
                ###################################################
                TL = np.zeros((0))
                for k in xrange(len(x)):
                    for l in xrange(int(round(y[k]*100))):
                        TL = np.append(TL, x[k])
                MTL_simulated[h, i, j] = TL.mean()
                
                ##############################################
                # calculate difference observed and simulated:
                ##############################################
                MTL_new = abs(MTL_data - TL.mean())
                if MTL_new<MTL_:
                    MTL[h] = TL.mean()
                    TLSTD[h] = TL.std()    
                MTL_ = MTL_new
        
        # record min and max simulated mean track lengh of each sample
        MTL_simulated_min[h] = MTL_simulated[h, :, :].min()
        MTL_simulated_max[h] = MTL_simulated[h, :, :].max()
    
        
        #################################################################################
        # calculate difference observed and simulated track lengths for each afta sample:
        #################################################################################
        if MTL_data >=  MTL_simulated[h, :, :].min() and MTL_data <=  MTL_simulated[h, :, :].max():
            MTL_diff[h] = 0   # i.e., observed data fall in range of simulated mean track lengths
        elif MTL_data<MTL_simulated[h, :, :].min():
            MTL_diff[h] = abs(MTL_data - MTL_simulated[h, :, :].min())  # observed mean length is lower than simulated length
        elif MTL_data>MTL_simulated[h, :, :].max():
            MTL_diff[h] = abs(MTL_data - MTL_simulated[h, :, :].max()) # observed mean length is higher than simulated length
        
    ##################################################################
    # calculate avg. difference mean and SD track lengths for all samples
    ######################################################################
    MTL_difference = 0
    TLSD_difference = 0
    for h in xrange(Nsamples):
        trackLengths_observed_ = trackLengths_observed[h]
        MTL_data = trackLengths_observed_.mean()
        TLSD_data = trackLengths_observed_.std()
        
        MTL_difference += MTL_diff[h]**2
        TLSD_difference += abs(TLSD_data - TLSTD[h])
    MTL_difference = math.sqrt(MTL_difference/Nsamples)   # rmse of mean difference
    TLSD_difference = TLSD_difference/Nsamples
    
    
    
    return MTL_difference, TLSD_difference, MTL_observed, MTL_simulated_min, MTL_simulated_max    
                         
    
def findCentralAges(AFTTable, sample_ids):
    
    print 'find central age data'
    
    AFT_central_age = np.zeros((len(sample_ids)))
    AFT_central_age_SE = np.zeros((len(sample_ids)))
    
    for i, sample in enumerate(sample_ids):
        index_ = np.where(AFTTable.cols.sample_id[:] == sample)[0]
        #condition = "sample_id  ==  '%s'" %sample
        
        AFT_central_age[i] = AFTTable.cols.central_age[index_]
        AFT_central_age_SE[i] = AFTTable.cols.central_age_SE[index_]
        
    return AFT_central_age, AFT_central_age_SE

def findBinomfitAges(binomFitTable, sample_ids):
    
    #print 'get binomfit ages from hdf5 database'
    
    AFT_Page = []
    AFT_Page_minConf = []
    AFT_Page_plusConf = []
    AFT_P_grains = []
    AFT_P_fraction = []
    
    for i, sample in enumerate(sample_ids):
        index_ = np.where(binomFitTable.cols.sample_id[:] == sample)[0]
        #condition = "sample_id  ==  '%s'" %sample

        AFT_Page_ = []
        AFT_Page_minConf_ = []
        AFT_Page_plusConf_ = []
        AFT_P_grains_ = []
        AFT_P_fraction_ = []
        
        if len(index_)>0:       
            # P1
            P1 = binomFitTable.cols.age1[index_]
            if P1> - 9999:
                AFT_Page_.append(P1)
                AFT_Page_minConf_.append(binomFitTable.cols.min_conf1[index_])
                AFT_Page_plusConf_.append(binomFitTable.cols.plus_conf1[index_])
                AFT_P_grains_.append(binomFitTable.cols.grains1[index_])
                AFT_P_fraction_.append(binomFitTable.cols.fraction1[index_])
            
            # P2
            P2 = binomFitTable.cols.age2[index_]
            if P2> - 9999:
                AFT_Page_.append(P2)
                AFT_Page_minConf_.append(binomFitTable.cols.min_conf2[index_])
                AFT_Page_plusConf_.append(binomFitTable.cols.plus_conf2[index_])
                AFT_P_grains_.append(binomFitTable.cols.grains2[index_])
                AFT_P_fraction_.append(binomFitTable.cols.fraction2[index_])
            
            # P3
            P3 = binomFitTable.cols.age3[index_]
            if P3> - 9999:
                AFT_Page_.append(P3)
                AFT_Page_minConf_.append(binomFitTable.cols.min_conf3[index_])
                AFT_Page_plusConf_.append(binomFitTable.cols.plus_conf3[index_])
                AFT_P_grains_.append(binomFitTable.cols.grains3[index_])
                AFT_P_fraction_.append(binomFitTable.cols.fraction3[index_])
        
        AFT_Page.append(AFT_Page_)
        AFT_Page_minConf.append(AFT_Page_minConf_)
        AFT_Page_plusConf.append(AFT_Page_plusConf_)
        AFT_P_grains.append(AFT_P_grains_)
        AFT_P_fraction.append(AFT_P_fraction_)
    
    return AFT_Page, AFT_Page_minConf, AFT_Page_plusConf, AFT_P_grains, AFT_P_fraction


def findAFTdata(AFTTable, well, removeCenozoicSamples=False,
                useGeotrackOnly=False):
    
    print 'read AFT age data from hdf5 database'
    
    ###################################
    # find AFT sample info (id, depth):
    ###################################
    condition = "well  ==  '%s'" %well
    condition_a = (AFTTable.cols.well[:] == well)
    condition_b = (AFTTable.cols.N_grains[:]>6)
    if useGeotrackOnly == True:
        condition_c = (AFTTable.cols.data_source[:] == 'Geotrack')
    else:
        condition_c = True
    sample_index = np.where(condition_a*condition_b*condition_c == True)[0]
    
    # remove Cenozoic samples:
    strat = AFTTable.cols.stratigraphy_code[:][sample_index]
    if removeCenozoicSamples == True:
        del_index = []
        for i, s in enumerate(strat):
            if s[0] == 'N': 
                del_index.append(i)
        if del_index !=  []:
            sample_index = np.delete(sample_index, del_index)
        
    # get AFT sample data:
    AFTsample_id = AFTTable.cols.sample_id[:][sample_index]
    AFTDepth_base = AFTTable.cols.depth_TVD_base[:][sample_index]
    AFTDepth_top = AFTTable.cols.depth_TVD_top[:][sample_index]
    Dpar_mean = AFTTable.cols.Dpar_mean[:][sample_index]
    Dpar_min = AFTTable.cols.Dpar_min[:][sample_index]
    Dpar_max = AFTTable.cols.Dpar_max[:][sample_index]
    
    # sort samples from shallow to deep:
    AFT_depth_index = np.argsort(AFTDepth_base)
    AFTsample_id_ = np.asarray(AFTsample_id)
    AFTsample_id_ = AFTsample_id_[AFT_depth_index]
    AFTsample_id = AFTsample_id_.tolist()
    AFTDepth_base = AFTDepth_base[AFT_depth_index]
    AFTDepth_top = AFTDepth_top[AFT_depth_index]
    Dpar_mean = Dpar_mean[AFT_depth_index]
    Dpar_min = Dpar_min[AFT_depth_index]
    Dpar_max = Dpar_max[AFT_depth_index]
    
    return  AFTsample_id, AFTDepth_base, AFTDepth_top,\
            Dpar_mean, Dpar_min, Dpar_max


def addPrePermErosion(subsidence_code, strat_units,
    cellSize_initial_all, cellSize_present_all, age_base, age_top,
    prePerm_erosion_start, prePerm_erosion_end, prePermEroded):
    
    
    print 'add late Carboniferous/early Permian erosion event'
    print prePermEroded
    
    ######################################################################################
    # add section that was eroded during the late - Carboniferous - early Permian unconformity
    # modifies the subsidence, stratigrphy code  lists and the cellsize and age arrays
    ######################################################################################
    
    if 'DCCU' in strat_units: topPermFm = 'DCCU'
    if 'DCDH' in strat_units: topPermFm = 'DCDH'
    if 'DCHL' in strat_units: topPermFm = 'DCHL'
    
    for i in xrange(len(strat_units)):
        if strat_units[i][:2] == 'DC': topPermFm = strat_units[i]
        
    prePermIndex = strat_units.index(topPermFm) + 1
    
    # update subsidence codes:
    subsidence_code.insert(prePermIndex, 'add_formation')
    subsidence_code.insert(prePermIndex + 1, 'erosion')
    
    # update stratigraphic codes table
    strat_units.insert(prePermIndex, (topPermFm + '-2'))
    strat_units.insert(prePermIndex + 1, ('-' + topPermFm + '-2'))
    
    # add eroded formations to initial cell size array:
    cellSize_initial_all = np.insert(cellSize_initial_all, prePermIndex, prePermEroded)
    cellSize_initial_all = np.insert(cellSize_initial_all, prePermIndex + 1, - prePermEroded)
    
    #
    cellSize_present_all = np.insert(cellSize_present_all, prePermIndex, prePermEroded)
    cellSize_present_all = np.insert(cellSize_present_all, prePermIndex + 1, - prePermEroded)
    
    # adjust ages of youngest Carboniferous formation:
    age_base = np.insert(age_base, prePermIndex, 306.5)
    age_top[prePermIndex - 1] = 306.5
    age_top = np.insert(age_top, prePermIndex, 299.)
    
    # age of (Stephanian) erosion phase:
    age_base = np.insert(age_base, prePermIndex + 1, prePerm_erosion_start)
    age_top = np.insert(age_top, prePermIndex + 1, prePerm_erosion_end)
    
    return subsidence_code, strat_units, cellSize_initial_all, cellSize_present_all, age_base, age_top


def addErosion(subsidence_code, strat_units, cellSize_initial_all,
    cellSize_present_all, age_base, age_top, erosion_stratCode,
    erosion_stratCode2, erosion_start, erosion_end, eroded_thickness):
    
    ######################################################################################
    # add eroded sections
    # modifies the subsidence, stratigrphy code  lists and the cellsize and age arrays
    ######################################################################################
    NerosionPhases = len(erosion_stratCode)
    
    for phaseNo in xrange(NerosionPhases):
        
        #print 'add erosion phase %s -  %s:' %(phaseNo + 1, erosion_stratCode[phaseNo])
        
        
        stratLn = len(erosion_stratCode[phaseNo])
        topFm = [] ; stratPresent = False
        
        # find eroded strat. member:
        for i in xrange(len(strat_units)):
            if strat_units[i][:stratLn] == erosion_stratCode[phaseNo]: topFm = strat_units[i] ; stratPresent = True
        
        # if strat. member not found, search for formation instead:
        if topFm == []:
            stratLn = len(erosion_stratCode2[phaseNo])
            for i in xrange(len(strat_units)):
                if strat_units[i][:stratLn] == erosion_stratCode2[phaseNo]: topFm = strat_units[i]
            
        if topFm !=  []:
                
            #print 'add fm:'
            
            if stratPresent == True:
                # insert new formation:
                    
                erosionIndex = strat_units.index(topFm) + 1

                # update subsidence codes:
                subsidence_code.insert(erosionIndex, 'add_formation')
                subsidence_code.insert(erosionIndex + 1, 'erosion')
                
                # update stratigraphic codes table
                strat_units.insert(erosionIndex, (topFm + '-2'))
                strat_units.insert(erosionIndex + 1, ('-' + topFm + '-2'))
                
                # add eroded formations to initial cell size array:
                cellSize_initial_all = np.insert(cellSize_initial_all, erosionIndex, eroded_thickness[phaseNo])
                cellSize_initial_all = np.insert(cellSize_initial_all, erosionIndex + 1, - eroded_thickness[phaseNo])
                
                #
                cellSize_present_all = np.insert(cellSize_present_all, erosionIndex, eroded_thickness[phaseNo])
                cellSize_present_all = np.insert(cellSize_present_all, erosionIndex + 1, - eroded_thickness[phaseNo])
                
                # add new ages of eroded section
                #start_S = age base_all[erosionIndex]
                midpoint = (age_base[erosionIndex - 1] + erosion_start[phaseNo])/2.
                
                age_top[erosionIndex - 1] = midpoint
                age_base = np.insert(age_base, erosionIndex, midpoint)    # timing of eroded section
                age_top = np.insert(age_top, erosionIndex, erosion_start[phaseNo]) # timing of eroded section  - > deposition up ot erosion start
                
                # add ages of erosion phase:
                age_base = np.insert(age_base, erosionIndex + 1, erosion_start[phaseNo]) # timing of erosion phase
                age_top = np.insert(age_top, erosionIndex + 1, erosion_end[phaseNo]) # timing of erosion phase

            else:
                
                # formation not present in well:
                
                erosionIndex = strat_units.index(topFm) + 1

                # update subsidence codes:
                subsidence_code.insert(erosionIndex, 'add_formation')
                subsidence_code.insert(erosionIndex + 1, 'erosion')
                
                # update stratigraphic codes table
                strat_units.insert(erosionIndex, (topFm + '-2'))
                strat_units.insert(erosionIndex + 1, ('-' + topFm + '-2'))
                
                # add eroded formations to initial cell size array:
                cellSize_initial_all = np.insert(cellSize_initial_all, erosionIndex, eroded_thickness[phaseNo])
                cellSize_initial_all = np.insert(cellSize_initial_all, erosionIndex + 1, - eroded_thickness[phaseNo])
                
                #
                cellSize_present_all = np.insert(cellSize_present_all, erosionIndex, eroded_thickness[phaseNo])
                cellSize_present_all = np.insert(cellSize_present_all, erosionIndex + 1, - eroded_thickness[phaseNo])
                
                # add new ages of eroded section
                midpoint = age_top[erosionIndex - 1] # start age of new unit equal to end age of underlying unit
                
                #age_top[erosionIndex - 1] = 
                age_base = np.insert(age_base, erosionIndex, midpoint)    # timing of eroded section
                age_top = np.insert(age_top, erosionIndex, erosion_start[phaseNo]) # timing of eroded section  - > deposition up ot erosion start
                
                # add ages of erosion phase:
                age_base = np.insert(age_base, erosionIndex + 1, erosion_start[phaseNo]) # timing of erosion phase
                age_top = np.insert(age_top, erosionIndex + 1, erosion_end[phaseNo]) # timing of erosion phase


                
    # 306.5 - 299 = age of stephanian
    #print 'stratigraphy after adding eroded fms'
    #print strat_units
    #print subsidence_code
    #print age_base
    #print age_top
    #exit()
    return subsidence_code, strat_units, cellSize_initial_all, cellSize_present_all, age_base, age_top



def gridRefinement(maxDuration_, maxCellSize, cellSize_initial_all,
    cellSize_present_all, strat_units, subsidence_code,
    age_base, age_top):
   
    print 'refine time steps'
   
    ############################################################################
    # subdivide tectonic time steps into smaller units if they exceed a certain thickness or duration
    ############################################################################
    Nsamples = len(cellSize_present_all)
    Ntectonic_steps = len(strat_units)
    i = 0 ; add = 0 ; erod = 0

    
    while i<Ntectonic_steps:
        
        duration = age_base[i] - age_top[i]
        
        maxDuration = maxDuration_
        
        if verbose:
            print i
            print strat_units[i]
            print cellSize_initial_all[i]
        
        if duration>maxDuration or cellSize_initial_all[i]>maxCellSize and strat_units[i][0] !=  '-' \
        and strat_units[i] !=  'FAULT' and strat_units[i][ - 2] !=  'f' and strat_units[i][ - 3] !=  'f':
            
            # calculate number of units to subdivide strat unit
            if (duration/maxDuration) >=  (cellSize_initial_all[i]/maxCellSize):
                divisions = int(round(duration/maxDuration) + 0.5)
            else:
                divisions = int(round(cellSize_initial_all[i]/maxCellSize))
            
            if divisions>1 and strat_units[i][0] !=  'F' and strat_units[i][0] !=  '-':
                
                duration_new = duration/divisions
                cellSize_new = cellSize_initial_all[i]/divisions
                cellSize_present_new = cellSize_present_all[i]/divisions
                strat_ = strat_units[i]
                subsidence_code_ = subsidence_code[i]
                #print 'subdividing unit %s - %s, %s in %s parts' %(i, strat_, subsidence_code[i], divisions)
                
                for j in xrange(0, divisions):
                    b = strat_ + '-f' + str(j)
                    if j == 0:
                        strat_units[i] = b
                        age_top[i] = age_base[i] - duration_new
                        cellSize_initial_all[i] = cellSize_new
                        cellSize_present_all[i] = cellSize_present_new
                        
                    else:
                        add += 1
                        strat_units.insert(i + 1, b)
                        subsidence_code.insert(i + 1, subsidence_code_)
                        cellSize_initial_all = np.insert(cellSize_initial_all, i + 1, cellSize_new)
                        cellSize_present_all = np.insert(cellSize_present_all, i + 1, cellSize_present_new)
                        age_base = np.insert(age_base, i + 1, age_base[i] - duration_new)
                        age_top = np.insert(age_top, i + 1, age_top[i] - duration_new)
                        i += 1
                        
                # check if unit is later eroded:
                strat_erosional = '-' + strat_
                if strat_erosional in strat_units[i:] and subsidence_code_ !=  'non_deposition':
                    k = strat_units.index(strat_erosional)
                    if subsidence_code[k] !=  'non_deposition':
                        #print 'subdividing corresponding erosional unit %s - %s, %s, in %s parts'\
                        #%(i, strat_erosional, subsidence_code[k], divisions)
                        # calculate new start and end time of deposition of unit:
                        duration = age_base[k] - age_top[k]
                        duration_new = duration/divisions
                        
                        # calculate new cell size:
                        cellSize_new = cellSize_initial_all[k]/divisions
                        cellSize_present_new = cellSize_present_all[k]/divisions
                        
                        # copy stratigraphy and subsidence code from unit:
                        strat_ = strat_units[k]
                        subsidence_code_ = subsidence_code[k]
                                    
                        # add new (subdivided) units:
                        for j in xrange(0, divisions):
                            b = strat_ + '-f' + str(j)
                            
                            if j == 0:
                                strat_units[k] = b
                                age_base[k] = age_top[k] + duration_new # start of new unit
                                age_base_ = age_base[k]   #
                                age_top_ = age_top[k]
                                cellSize_initial_all[k] = cellSize_new    # new initial cell size:
                                cellSize_present_all[k] = cellSize_present_new    # new present day cell size:
                            else:
                                erod += 1
                                strat_units.insert(k, b)
                                subsidence_code.insert(k, subsidence_code_)
                                cellSize_initial_all = np.insert(cellSize_initial_all, k, cellSize_new)
                                cellSize_present_all = np.insert(cellSize_present_all, k, cellSize_present_new)
                                age_base = np.insert(age_base, k, age_base_ + j*duration_new)
                                age_top = np.insert(age_top, k, age_top_ + j*duration_new)
                                
        i += 1    
        Ntectonic_steps = len(subsidence_code)
        
    return cellSize_initial_all, cellSize_present_all, strat_units, subsidence_code, age_base, age_top


def calculateDecompaction(cellSize_start, fractionArray, Phi0, Phic,
    NcompactionIterations):
    
    # calculate depths of cell midpoints    
    cellDepths_start = calculateDepths_mid(cellSize_start)
    # calculate porosity
    phiArray_start = porosityDepthFunc(cellDepths_start, fractionArray, Phi0, Phic)
    # start with decompacted cell size = present cell size
    cellSize_final = cellSize_start.copy()    
    
    for j in xrange(NcompactionIterations):
        
        # calculate porosity after removing overlying units:
        phiArray_final = porosityDepthFunc(cellSize_final/2.0, fractionArray, Phi0, Phic)
        # calculate cell size from comparison surface porosity and present day porosity:
        cellSize_initial = calculateCompactedThickness(cellSize_start, phiArray_start, phiArray_final)
    
    return cellSize_initial

def checkBurial(strat_units, subsidence_code, cellSize_present,
    cellSize_initial):
    
    for i in xrange(len(strat_units)):
        print '%s: %s, %s, init th =  %.0f, present th = %.0f' %(i, strat_units[i], subsidence_code[i], cellSize_initial[i], cellSize_present[i])
    
    return

def prepareWellData(
        h5file, well, maxDuration, maxCellSize, 
        lateCret_erosion_start, lateCret_erosion_end, 
        lateCret_erosion, erosion_stratCode, erosion_stratCode2,
        erosion_start, erosion_end, eroded_thickness, 
        prePerm_erosion_start, prePerm_erosion_end, prePermEroded, 
        Phi0, Phic, maxATthickness_ ):

        
    print 'prepare 1D subsidence history data'
    
    #################################
    # find stratigraphy of wells:
    #################################
    
    # find stratigraphic ages of well using lookup table in hdf5 database
    wellStratTable = h5file.root.strat.strat_mod
    #stratCodesTable = h5file.root.strat.strat_petroprob
    stratCodesTable = h5file.root.strat.formation_descriptions

    # find well stratigraphy, ages and lithology in hdf5 database:
    # new alternative:
    stratCodesTable = h5file.root.strat.formation_descriptions
    
    wellStratCode, wellStratDepth_baseTVD, wellStratDepth_topTVD =\
        findStrat(wellStratTable, well)
    age_base_well, age_top_well = findStratAges(
                                    wellStratCode, stratCodesTable)
    
    stratCodesTable = h5file.root.strat.formation_descriptions
    fractionArray = findLithologies_1D(wellStratCode, stratCodesTable)

    # set grid cell size arrays from present day formation depths
    cellSize_present_well = setGridCells(wellStratDepth_baseTVD[:len(wellStratCode)])
    
    # find initial cell size, using simple decompaction
    print 'decompaction'
    NcompactionIterations = 5
    cellSize_initial_well = calculateDecompaction(
                cellSize_present_well.copy(), fractionArray, 
                Phi0, Phic, NcompactionIterations)

    # print decompaction results
    if verbose == True:
        for i in xrange(len(cellSize_initial_well)):
        
            print 'decomp -  %s, size: %0.0f to %0.0f'\
            %(wellStratCode[i], cellSize_present_well[i], cellSize_initial_well[i])
        
    # set up time steps, using stratigraphy (subsidence or non_deposition)  + info on erosion phases
    print 'arrange tectonic time steps:'

    # create list with tectonic labels for each time step ( - > 'subsidence', 'erosoion' or 'non_deposition'):
    Ncells = len(cellSize_initial_well)
    subsidence_code = ['subsidence']*Ncells   # initial: subsidence for each stratigrpahic unit present
    strat_units = wellStratCode[:: - 1]           # copy well stratigraphy
    age_base = age_base_well[:: - 1] ; age_top = age_top_well[:: - 1] # copy stratigraphic ages
    cellSize_initial = cellSize_initial_well[:: - 1] ; cellSize_present = cellSize_present_well[:: - 1]   # copy cell sizes

    prePerm = True ; totalATerosion = 0 ; totalSLerosion = 0

    ########################################################################
    # add eroded formations (Jurassic Altena group +  lower Cretaceous SLDNA)
    ########################################################################
    erodedATTable = h5file.root.subsidence_model.eroded_thickness_Altena
    memberTable = h5file.root.subsidence_model.member_thickness_Altena    
        
    # add eroded Altena Gp. (AT) formations, only if lower Cretaceous SLDNA also eroded 
    print 'add eroded formations'
    if 'SLDNA' not in strat_units:
        maxErosion = lateCret_erosion
        if maxErosion>maxATthickness_:
            maxErosion = maxATthickness_
        subsidence_code, strat_units, cellSize_initial, cellSize_present, age_base, age_top\
        , lateCret_erosion_base, lateCret_erosion_top, totalATerosion = \
        addErodedAltenaFms(well, erodedATTable, memberTable, subsidence_code, strat_units\
        , cellSize_initial, cellSize_present, age_base, age_top, lateCret_erosion_start, lateCret_erosion_end\
        , maxErosion)
    else:
        lateCret_erosion_base = 0 ; lateCret_erosion_top = 0 ; totalATerosion = 0

    print '\tEroded thickness Altena fm: %0.0f' %totalATerosion
        
    subsidence_code, strat_units, cellSize_initial, cellSize_present, age_base, age_top, totalSLerosion = \
    addErodedSLDNAfm(subsidence_code, strat_units, cellSize_initial.copy(), cellSize_present.copy()\
    , age_base, age_top, lateCret_erosion_start, lateCret_erosion_end, lateCret_erosion_base\
    , lateCret_erosion_top, lateCret_erosion, totalATerosion)

    print '\tEroded thickness Schieland fm %0.0f ' %totalSLerosion
    
    # add pre-Permian erosion event, if Carboniferous fromation penetrated by well:
    if 'DCCU' in strat_units or 'DCDH' in strat_units or 'DCHS' in strat_units:
        print 'add permian erosion event'
        subsidence_code, strat_units, cellSize_initial,\
        cellSize_present, age_base, age_top = addPrePermErosion(
                subsidence_code, strat_units,
                cellSize_initial, cellSize_present, age_base, age_top,
                prePerm_erosion_start, prePerm_erosion_end,
                prePermEroded)
        prePerm = True
    else:
        prePerm = False
    
    
    
    ###################################      
    # add additional erosion events
    ###################################
    print 'add additional erosion events'
    subsidence_code, strat_units, cellSize_initial, cellSize_present, age_base, age_top = \
    addErosion(subsidence_code, strat_units, cellSize_initial.copy(), cellSize_present.copy(), age_base, age_top\
    , erosion_stratCode, erosion_stratCode2, erosion_start, erosion_end, eroded_thickness)
            
    Ntectonic_steps = len(subsidence_code)

    ########################################
    # add non - erosional unconformities to subsidence code list:
    ########################################
    print 'add non - deposition phases'
    subsidence_code, strat_units, cellSize_initial, cellSize_present, age_top, age_base\
     = add_nonDeposition(subsidence_code, strat_units, cellSize_initial.copy(), cellSize_present.copy(), age_top, age_base)
    Ntectonic_steps = len(subsidence_code)

    # check to correct eroneous tectonic time step ages
    age_base, age_top = correctAges(age_base.copy(), age_top.copy())

    
    ##########################
    # subdivide steps:
    ##########################
    cellSize_initial, cellSize_present, strat_units, subsidence_code, age_base, age_top = \
    gridRefinement(maxDuration, maxCellSize, cellSize_initial.copy(), cellSize_present.copy()\
    , strat_units, subsidence_code, age_base, age_top)

    if verbose == True:
        print 'check burial sequence:'
        showCellSizes_present(cellSize_present, cellSize_initial, strat_units)
    
    
    ###########################################
    # test decompaction after grid refinement:
    ###########################################
    print 'adjust decompaction for erosional phases'
    useSecondDecomp = True
    if useSecondDecomp == True:
        
        # leave out fms that will be eroded later:
        non_erosion_index = np.zeros((0), dtype = int)
        for i in xrange(len(subsidence_code) - 1, - 1, - 1):
            if subsidence_code[i] !=  'add_formation' and subsidence_code[i] !=  'erosion' and subsidence_code[i] !=  'non_deposition':
                non_erosion_index = np.append(non_erosion_index, i)
        
        # select present and initial cell sizes of non - eroded fm's
        cellSize_present_nonEroded = cellSize_present[non_erosion_index]
        cellSize_initial_nonEroded = cellSize_initial[non_erosion_index]
        
        # find non - eroded strat units:      
        strat_units_nonEroded = [] ; strat_units_clean = []
        for j, s in enumerate(strat_units):
            if j in non_erosion_index:
                strat_units_nonEroded.append(s)  
                strat_units_clean.append(s)      
            else:
                strat_units_clean.append('NUCT')      
        strat_units_nonEroded = strat_units_nonEroded[:: - 1]
        
        # find lithologies of non - eroded fms
        fractionArray_nonEroded = findLithologies_1D(strat_units_nonEroded, stratCodesTable)
    
        print '\tdecompaction'
        cellSize_initial_nonEroded = calculateDecompaction(cellSize_present_nonEroded.copy(), fractionArray_nonEroded\
        , Phi0, Phic, NcompactionIterations)
        
        cellSize_initial[non_erosion_index] = cellSize_initial_nonEroded.copy()
        

        ###################################
        # correct eroded section thickness:
        ###################################
        # find pre-inversion stratigraphy:       
        strat_eroded = [] ; eroded_index = []
        cellSize_eroded_lst = [] ; cellSize_initial_eroded_lst = []
        count = 0
        
        # skip carboniferous/early Permian
        if 'DC' in strat_units[count]:
            while 'DC' in strat_units[count]: count += 1
        while strat_units[count][0] !=  '-' and strat_units[count][1] !=  'D':
            if subsidence_code[count] == 'add_formation':
                strat_eroded.append(strat_units[count])
                cellSize_eroded_lst.append(cellSize_present[count])
                cellSize_initial_eroded_lst.append(cellSize_initial[count])
                eroded_index.append(count)
            count += 1
        cellSize_present_eroded = np.asarray(cellSize_eroded_lst)[:: - 1]
        cellSize_initial_eroded = np.asarray(cellSize_initial_eroded_lst)[:: - 1]
        strat_eroded = strat_eroded[:: - 1]
        
        # calculate lithologies:
        fractionArray_eroded = findLithologies_1D(strat_eroded, stratCodesTable)
        
        print '\tdecompacting eroded section'
        cellSize_initial_eroded = calculateDecompaction(cellSize_present_eroded.copy(), fractionArray_eroded\
        , Phi0, Phic, NcompactionIterations)

        print '\tsum eroded section, present = %0.0f, initial = %0.0f' %(cellSize_present_eroded.sum(), cellSize_initial_eroded.sum())
        
        # update eroded section thickness:
        cellSize_initial[eroded_index] = cellSize_initial_eroded
        
        ############################################
        # test deeper burial before erosion events
        ############################################
        # find pre-inversion stratigraphy:       
        strat_preInversion = [] ; cellSize_preInversion_lst = [] ; cellSize_initial_preInversion_lst = [] ; pre_inversion_index = []
        count = 0
        # skip pre perm erosion:
        if 'DC' in strat_units[count]:
            while 'DC' in strat_units[count]: count += 1
        # find pre-inversion strat succession:
        while strat_units[count][0] !=  '-' or 'DC' in strat_units[count]:
            if subsidence_code[count] !=  'non_deposition':
                strat_preInversion.append(strat_units[count])
                cellSize_preInversion_lst.append(cellSize_present[count])
                cellSize_initial_preInversion_lst.append(cellSize_initial[count])
                pre_inversion_index.append(count)
            count += 1
        #print 'pre inversion strat =  %s' %(strat_preInversion)
        strat_preInversion = strat_preInversion[:: - 1]
        cellSize_present_preInversion = np.asarray(cellSize_preInversion_lst)[:: - 1]
        cellSize_initial_preInversion = np.asarray(cellSize_initial_preInversion_lst)[:: - 1]
        
        # calculate lithologies:
        fractionArray_preInversion = findLithologies_1D(strat_preInversion, stratCodesTable)
        
        print 'decompacting pre-erosion stratigraphic units'
        cellSize_initial_preInversion = calculateDecompaction(cellSize_present_preInversion.copy(), fractionArray_preInversion\
        , Phi0, Phic, NcompactionIterations)

        # compare pre-inversion and final calibrated initial cell sizes:
        print 'updating strat units for pre-inversion deep burial:'
        for i in xrange(len(strat_preInversion)):
            if strat_preInversion[i] in strat_units:
                index_ = strat_units.index(strat_preInversion[i])
                # update initial cell sizes
                if cellSize_initial_preInversion[i]>cellSize_initial[index_]:
                    cellSize_initial[index_] = cellSize_initial_preInversion[i]
                    pass
                if verbose == True:
                    print '%s - %s, pre-inv cellsize = %0.0f, post - inv cellsize = %0.0f' %(strat_preInversion[i], strat_units[index_], cellSize_initial_preInversion[i]\
                    , cellSize_initial[index_])
        
    fractionArray_allSteps = []
  
    return cellSize_initial, cellSize_present, strat_units\
    , subsidence_code, age_base, age_top\
    , stratCodesTable, wellStratCode, wellStratDepth_topTVD, wellStratDepth_baseTVD, prePerm\
    , totalATerosion, totalSLerosion


def simulateStratUnits(strat_units, subsidence_code):
    
    Nstrat = len(strat_units)
    Ntectonic_steps = len(subsidence_code)
    #print strat_units[0]
    #exit()
    strat_steps = [[strat_units[0]]]
    
    for i in xrange(1, Ntectonic_steps):
        strat_steps_local = strat_steps[i - 1][:]
        if subsidence_code[i] == 'subsidence' or subsidence_code[i] == 'add_formation':                
            strat_steps_local.insert(0, strat_units[i])
        elif subsidence_code[i] == 'erosion':
            if strat_units[i][0] == '-':
                # find eroded section
                if '-' in strat_units[i]:
                    s = strat_units[i][1:]
                else:
                    print 'error, no - in subsidence code'
                    print '%s , %s' %(strat_units[i], subsidence_code[i])
                    exit()
                if s in strat_steps_local[0]:
                    del(strat_steps_local[0])
                else:
                    print 'error erosion algorithm, eroded unit %s not top unit %s' %(s, strat_steps[i][0])
                    exit()
        
        strat_steps.append(strat_steps_local)
        
    return strat_steps


def simulateBurial(burialParam, subsidence_code, Nstrat):
    
    #cellSize_steps = simulateBurial(cellSize_initial, subsidence_code)
    #cellSize_initial_steps = simulateBurial(cellSize_initial, subsidence_code)
    
    
    #Nstrat = len(burialParam)
    Ntectonic_steps = len(subsidence_code)
    
    burialParam_steps = np.zeros((Ntectonic_steps, Nstrat))
    
    for i in xrange(Ntectonic_steps):
    
        if subsidence_code[i] == 'subsidence' or subsidence_code[i] == 'add_formation':                
            # insert new stratigraphic unit and surface temperature
            burialParam_steps[i][1:] = burialParam_steps[i - 1][: - 1]
            try:
                burialParam_steps[i][0] = burialParam[i]
            except:
                print len(burialParam)
                print i
                exit()
        elif subsidence_code[i] == 'erosion':
            burialParam_steps[i][: - 1] = burialParam_steps[i - 1][1:]
        else:
            burialParam_steps[i] = burialParam_steps[i - 1]
            
    return burialParam_steps


def findLithologies(Nstrat, strat_steps, stratCodesTable):
    #fractionArray = findLithologies(stratCodes, stratCodesTable)
    
    ## find formation fractions of each stratigraphic interval in a well (wellStratTable) from hdf5 database of strat. units and lithology 
    Ntectonic_steps = len(strat_steps)
    Nlithologies = 5
    fractionArray = np.zeros((Ntectonic_steps, Nstrat, Nlithologies))
    
    sand = stratCodesTable.cols.sand[:]
    silt = stratCodesTable.cols.silt[:]
    congl = stratCodesTable.cols.congl[:]
    clay = stratCodesTable.cols.clay[:]
    shale = stratCodesTable.cols.shale[:]
    lime = stratCodesTable.cols.lime[:]
    carbonate = stratCodesTable.cols.carbonate[:]
    chalk = stratCodesTable.cols.chalk[:]
    oolite = stratCodesTable.cols.oolite[:]
    dolomite = stratCodesTable.cols.dolomite[:]
    marl = stratCodesTable.cols.marl[:]
    coal = stratCodesTable.cols.coal[:]
    organic = stratCodesTable.cols.organic[:]
    halite = stratCodesTable.cols.halite[:]
    anhydrite = stratCodesTable.cols.anhydrite[:]
    
    strat = stratCodesTable.cols.code_fm[:]
    
    #stratIndex = np.where(stratCodesTable.cols.code_fm[:] == member_code)[0][0]
    for i in xrange(Ntectonic_steps):
            
        for counter, member_code in enumerate(strat_steps[i]):
            if member_code !=  []:
                if member_code[0] == '-' :
                    print 'warning: eroded member %s in strat column' %member_code 
                    member_code = member_code[1:]
            
                for j in xrange(2, len(member_code) + 1):
                    if member_code[:j] in strat[:]:
                        Nchar = j
                try:
                    member_code = member_code[:Nchar]
                except:
                    print 'error, strat unit not in database'
                    print member_code
                    print strat_steps[i]
                    print Nchar
                    print stratCodesTable.cols.code_fm[:]
                    exit()
                
                stratIndex = np.where(strat == member_code)[0]
                
                fractionArray[i, counter, 0] = sand[stratIndex]
                fractionArray[i, counter, 0] += silt[stratIndex]
                fractionArray[i, counter, 0] += congl[stratIndex]
                
                fractionArray[i, counter, 1] = clay[stratIndex]
                fractionArray[i, counter, 1] += shale[stratIndex]
                
                fractionArray[i, counter, 2] = lime[stratIndex]
                fractionArray[i, counter, 2] += carbonate[stratIndex]
                fractionArray[i, counter, 2] += oolite[stratIndex]
                fractionArray[i, counter, 2] += chalk[stratIndex]
                fractionArray[i, counter, 2] += dolomite[stratIndex]
                fractionArray[i, counter, 2] += marl[stratIndex]
                
                fractionArray[i, counter, 3] = coal[stratIndex]
                fractionArray[i, counter, 3] += organic[stratIndex]
                
                fractionArray[i, counter, 4] = halite[stratIndex]
                fractionArray[i, counter, 4] += anhydrite[stratIndex]
            
                # convert percentages to fractions
                fractionArray[i, counter] = fractionArray[i, counter]/100.0
                                
                # correct for 0 or >1 formation fractions
                if fractionArray[i, counter, :].sum()>1.:
                    print 'warning, >100% lithology info for unit "%s"  ' %member_code
                    fractionArray[i, counter, :] = fractionArray[i, counter, :]/(fractionArray[i, counter, :].sum())
                if fractionArray[i, counter, :].sum() == 0:
                    fractionArray[i, counter, :4] = np.array([0.5, 0.5, 0., 0.])
                    print 'warning, no lithology info for unit "%s"  ' %member_code
        
    return fractionArray


def findLithologies_1D(strat_steps, stratCodesTable):
    #fractionArray = findLithologies(stratCodes, stratCodesTable)
    
    ## find formation fractions of each stratigraphic interval in a well (wellStratTable) from hdf5 database of strat. units and lithology 
    Nstrat = len(strat_steps)
    Nlithologies = 5
    fractionArray = np.zeros((Nstrat, Nlithologies))
    
    sand = stratCodesTable.cols.sand[:]
    silt = stratCodesTable.cols.silt[:]
    congl = stratCodesTable.cols.congl[:]
    clay = stratCodesTable.cols.clay[:]
    shale = stratCodesTable.cols.shale[:]
    #blackclay = stratCodesTable.cols.blackclay[:]
    lime = stratCodesTable.cols.lime[:]
    carbonate = stratCodesTable.cols.carbonate[:]
    chalk = stratCodesTable.cols.chalk[:]
    oolite = stratCodesTable.cols.oolite[:]
    dolomite = stratCodesTable.cols.dolomite[:]
    marl = stratCodesTable.cols.marl[:]
    coal = stratCodesTable.cols.coal[:]
    organic = stratCodesTable.cols.organic[:]
    halite = stratCodesTable.cols.halite[:]
    anhydrite = stratCodesTable.cols.anhydrite[:]
    
    strat = stratCodesTable.cols.code_fm[:]
        
    for counter, member_code in enumerate(strat_steps):
        if member_code !=  []:
            if member_code[0] == '-' :
                print 'warning: eroded member %s in strat column' %member_code 
                member_code = member_code[1:]
        
            for j in xrange(2, len(member_code) + 1):
                if member_code[:j] in strat[:]: Nchar = j 
            try:
                member_code = member_code[:Nchar]
            except:
                print 'error, strat unit not in database'
                print member_code
                print strat_steps[i]
                print Nchar
                print stratCodesTable.cols.code_fm[:]
                exit()
            
            stratIndex = np.where(strat == member_code)[0]
            
            try:
                fractionArray[counter, 0] = sand[stratIndex]
            except:
                print 'error, unit %s not found in lithology database' %(member_code)
                exit()
                
            fractionArray[counter, 0] += silt[stratIndex]
            fractionArray[counter, 0] += congl[stratIndex]
            
            fractionArray[counter, 1] = clay[stratIndex]
            fractionArray[counter, 1] += shale[stratIndex]
            
            fractionArray[counter, 2] = lime[stratIndex]
            fractionArray[counter, 2] += carbonate[stratIndex]
            fractionArray[counter, 2] += oolite[stratIndex]
            fractionArray[counter, 2] += chalk[stratIndex]
            fractionArray[counter, 2] += dolomite[stratIndex]
            fractionArray[counter, 2] += marl[stratIndex]
            
            fractionArray[counter, 3] = coal[stratIndex]
            fractionArray[counter, 3] += organic[stratIndex]
            
            fractionArray[counter, 4] = halite[stratIndex]
            fractionArray[counter, 4] += anhydrite[stratIndex]
        
            # convert percentages to fractions
            fractionArray[counter] = fractionArray[counter]/100.0
            
            # correct for 0 or >1 formation fractions
            if fractionArray[counter, :].sum()>1.:
                print 'warning, >100% lithology info for unit "%s"  ' %member_code
                fractionArray[counter, :] = fractionArray[counter, :]/(fractionArray[counter, :].sum())
            if fractionArray[counter, :].sum() == 0:
                fractionArray[counter, :4] = np.array([0.5, 0.5, 0., 0.])
                print 'warning, no lithology info for unit "%s"  ' %member_code
    
    return fractionArray


def porosityDepthFunc(depth, fractionArray, Phi0, Phic):
    
    """
    calculation of porosity using an exponential porosity - depth function
    """
    
    # calculate porosities
    try:
        Ntectonic_steps, Nstrat = depth.shape
        transformArray = False
    except:
        depth.shape = 1, - 1
        a, b = fractionArray.shape
        fractionArray.shape = 1, a, b
        Ntectonic_steps, Nstrat = depth.shape
        transformArray = True
        
    #porosity = np.zeros((Ntectonic_steps, Nstrat))
                
    sand = fractionArray[:, :, 0]
    shale = fractionArray[:, :, 1]
    carbonate = fractionArray[:, :, 2]
    coal = fractionArray[:, :, 3]
    evaporite = fractionArray[:, :, 4]
            
    # calculate porosity
    porosity = Phi0[0]*np.exp( - Phic[0]*depth)*shale\
     + Phi0[1]*np.exp( - Phic[1]*depth)*sand\
     + Phi0[2]*np.exp( - Phic[2]*depth)*carbonate\
     + Phi0[3]*coal\
     + Phi0[4]*evaporite           
    
    if transformArray:
        # transform back to 1D array
        porosity.shape =  - 1
        depth.shape =  - 1
        a, b, c = fractionArray.shape
        fractionArray.shape = b, c
        
    return porosity


def calculateDepths(cellSize):
    
    # calculate cell depths from a 2D array of cell dimensions
    Ntectonic_steps, Ncells = cellSize.shape
    cellDepths_base = np.zeros((Ntectonic_steps, Ncells))
    cellDepths_top = np.zeros((Ntectonic_steps, Ncells))
    cellDepths_top[:, 0] = 0
    for i in xrange(Ntectonic_steps):
        for j in range(Ncells):
            cellDepths_base[i, j] = cellSize[i, :(j + 1)].sum()
            if j>1:
                cellDepths_top[i, j] = cellDepths_base[i, j - 1]
    return cellDepths_base, cellDepths_top


def calculateDepths_mid(cellSize):
    
    # calculate cell depth midpoints from an 1D array of cell dimensions
    try:
        Ntectonic_steps, Ncells = cellSize.shape
        transformArray = False
    except:
        # transform 1D array to 2D:
        cellSize.shape = 1, - 1
        Ntectonic_steps, Ncells = cellSize.shape
        transformArray = True
    
    cellDepths_mid = np.zeros((Ntectonic_steps, Ncells))
    cellDepths_mid[:, 0] = cellSize[:, 0]/2.
    #for i in xrange(Ntectonic_steps):
    #    for j in range(1, Ncells):
    #            cellDepths_mid[i, j] = cellSize[i, :j].sum() + (cellSize[i, j]/2.)
    #cellDepths_mid2 = cellDepths_mid.copy()
    for j in range(1, Ncells):
        cellDepths_mid[:, j] =  np.sum(cellSize[:, :j],axis=1) +\
                                (cellSize[:, j]/2.)
    
    if transformArray:
        # transform arrays back to 1D
        #print cellDepths_mid
        cellSize.shape =  - 1
        cellDepths_mid.shape =  - 1
        
    return cellDepths_mid


def calculateDepths_base(cellSize):
    
    # calculate cell depth midpoints from an 1D array of cell dimensions
    try:
        Ntectonic_steps, Ncells = cellSize.shape
        transformArray = False
    except:
        # transform 1D array to 2D:
        cellSize.shape = 1, - 1
        Ntectonic_steps, Ncells = cellSize.shape
        transformArray = True
    
    cellDepths_base = np.zeros((Ntectonic_steps, Ncells))
    for i in xrange(Ntectonic_steps):
        for j in range(0, Ncells):
                #cellDepths[j] = cellDepths[j-1] + (cellSize[j-1]/2 + cellSize[j]/2)
                cellDepths_base[i, j] = cellSize[i, :j].sum() + (cellSize[i, j])
    
    if transformArray:
        # transform arrays back to 1D
        #print cellDepths_mid
        cellSize.shape =  - 1
        cellDepths_base.shape =  - 1
    return cellDepths_base


def calculateCompactedThickness(thickness_old, phi_old, phi_new):
    
    thickness_new = thickness_old/(( - phi_new + 1)/( - phi_old + 1)) #
            
    return thickness_new


def simulateCompaction(cellSize_initial, fractionArray, Phi0, Phic):

    #porosity, cellSize = simulateCompaction(cellSize_initial, fractionArray\
    #, Phi0, Phic, NcompactionIterations)

    Ntectonic_steps, Nstrat = cellSize_initial.shape
    #porosity = np.ones((Ntectonic_steps, Nstrat))*0.5
    
    # calculate initial estimate of new cell depths:
    cellDepths_compact = calculateDepths_mid(cellSize_initial)
    
    # calculate initial porosity (i.e. strat unit  = at surface)
    porosity_initial\
     = porosityDepthFunc(cellSize_initial/2.0, fractionArray, Phi0, Phic)
    
    porosity = porosity_initial.copy()
    
    # iterate compacted cell size:
    NcompactionIterations = 5
    for j in xrange(NcompactionIterations):
    
        # new porosity:
        porosity_compacted\
         = porosityDepthFunc(cellDepths_compact, fractionArray,
                                Phi0, Phic)
        
        # make sure porosity does not increase
        #(i.e., no change in porosity during erosion and reburial):
        porosity_compacted[np.where(porosity_compacted>porosity)] =\
                    porosity[np.where(porosity_compacted>porosity)]
        
        # calculate compacted cell size
        cellSize_compacted = calculateCompactedThickness(
                                cellSize_initial, porosity_initial,
                                porosity_compacted )
        
        # get cell depths from cell sizes:
        cellDepths_compact = calculateDepths_mid(cellSize_compacted)
    
    cellSize = cellSize_compacted.copy()
    porosity = porosity_compacted.copy()
    
    return porosity, cellSize
    

def calculateThermalConductivity(fractionArray, porosity,
    Kclay, Ksand, Kchalk, Korg, Kevap, Kwater):
    #Kbulk = calculateThermalConductivity(fractionArray, porosity, Kclay, Ksand, Kchalk, Korg, Kevap, Kwater)
    
    # calculate bulk conductivity vector from lithology dependent matrix conductivity values 
    
    # KBulk = array to be used for results (depth top, depth bottom)
    # fractionArray =  fractions of lithology (shale, sand, limestone, organic matter/lignite/coal, evaporite)
    # Ksand, Kchalk, Korg, Kwater = thermal conductivities
    
    Ntectonic_steps, Nstrat, Nlithologies = fractionArray.shape
    Kbulk = np.zeros((Ntectonic_steps, Nstrat))
    
    ############          
    for j in xrange(Ntectonic_steps):
        for i in xrange(Nstrat):
            
            sand, shale, carbonate, coal, evaporite = fractionArray[j, i]
            
            # calculate matrix thermal conductivity from lithology - dependent conductivity values and lithology - depth fractions file
            Kmatrix = Kclay*shale + Ksand*sand + Kchalk*carbonate + Korg*coal + Kevap*evaporite 

            # calculate bulk thermal conductivity form porosity and matrix and water thermal conductivity
            Kbulk[j, i] = (Kwater**porosity[j, i])*(Kmatrix**(1 - porosity[j, i]))
    return Kbulk
    

def setSurfaceTemp(surfaceTempHist, surfaceTempHistMarkers,
                    strat_units):
    
    Ntectonic_steps = len(strat_units)
    surfaceTemp = np.zeros((Ntectonic_steps))
    for i in xrange(Ntectonic_steps):
        s = strat_units[i]
        if s[0] == '-': s = s[1:3]
        else: s = s[:2]
        if s in surfaceTempHistMarkers:
            surfaceTemp[i] = surfaceTempHist[surfaceTempHistMarkers.index(s)]
        else:
            print 'error, unit %s not in surface temp array' %s
            print surfaceTempHistMarkers
            exit()
    return surfaceTemp


def setHeatFlow(heatflowMarkers, heatflowHist, strat_units):
    #heatFlow = setHeatFlow(heatflowMarkers, heatflowHist, strat_units)
    # get bottom heatflow        
    Ntectonic_steps = len(strat_units)
    heatFlow = np.zeros((Ntectonic_steps))
    for i in xrange(Ntectonic_steps):
        for j in xrange(len(heatflowMarkers)):
            if heatflowMarkers[j] == 'all':
                heatFlow[i] = heatflowHist[j]/1000.
            else:
                a = len(heatflowMarkers[j])
                if strat_units[i][:a] == heatflowMarkers[j]:
                    heatFlow[i] = heatflowHist[j]/1000
    
    return heatFlow

def changeCellThickness(cellSize_initial_steps, strat_units,
    strat_steps, changeThickness_timing, changeThickness_units,
    changeThickness_value):
#cellSize_initial_steps = changeCellThickness(cellSize_initial_steps, strat_units, strat_steps\
#, changeThickness_timing, changeThickness_units, changeThickness_value)

    Ntectonic_steps, Nstrat = cellSize_initial_steps.shape
    
    # find time steps for thickness changes:
    changeIndex_timing = []
    Nsteps = len(changeThickness_timing)
    changeIndex_timing2 = np.ones((Nsteps), dtype = int)*99999
    faultMovement = False
    for k in xrange(Ntectonic_steps):
        for i in xrange(Nsteps):
            if changeThickness_timing[i] in strat_units[k] and strat_units[k][0] !=  '-':
                if i == 0 and changeIndex_timing2[i] == 99999:
                    changeIndex_timing2[0] = k
                else:
                    changeIndex_timing2[i] = k
                faultMovement = True
     #       else:
     #           print 'unit %s not in %s' %(changeThickness_timing[i], strat_units[k])
    
    #exit()
                
    if faultMovement:
        # calculate thickness change for each timestep
        for i in xrange(1, Nsteps):
            Nsteps_dummy = changeIndex_timing2[i] - changeIndex_timing2[i - 1]
            if i == 1:
                thicknessChange_distributed = np.linspace(changeThickness_value[i - 1], changeThickness_value[i], Nsteps_dummy)
            else:
                thicknessChange_distributed = np.concatenate((thicknessChange_distributed\
                , np.linspace(changeThickness_value[i - 1], changeThickness_value[i], Nsteps_dummy)))
        Nsteps = changeIndex_timing2[-1] - changeIndex_timing2[0]
        
        print 'fault movement distributed over %i tectonic steps:' %Nsteps
        print strat_units[changeIndex_timing2[0]:changeIndex_timing2[-1]]
        #print strat_units[changeIndex_mid]
        #print strat_units[changeIndex_end]
        
        # find strat. units to change
        for k in xrange(changeIndex_timing2[0], changeIndex_timing2[-1]):
            changeIndex = []
            # check if strat unit is present:
            for l in xrange(len(strat_steps[k])):
                if changeThickness_units in strat_steps[k][l]:
                    changeIndex.append(l)
                    #print 'found unit %s in strat unit %s' %(changeThickness_units, strat_steps[k][l])
            # change thickness of units:
            if changeIndex !=  []:
                #exit()
                # distribute thickness change over affected units:
                Nthickness_changes = len(changeIndex)
                changeThickness_value_adj = thicknessChange_distributed[k - changeIndex_timing2[0]]/Nthickness_changes
                changeThickness = cellSize_initial_steps[k, changeIndex] + changeThickness_value_adj
                changeThickness[np.where(changeThickness<1)] = 1
                
                if verbose:
                    print 'strat step: %s' %strat_units[k]
                    print '- - - old thickness: ', cellSize_initial_steps[k, changeIndex].sum()
                
                cellSize_initial_steps[k, changeIndex] = changeThickness
                
                if verbose:
                    print '- - - changed cell size of units below with %s' %changeThickness_value_adj
                    print '- - - ', np.asarray(strat_steps[k])[changeIndex]
                    print '- - - new thickness: ', cellSize_initial_steps[k, changeIndex].sum()
    
    return cellSize_initial_steps


def setTimeSteps(age_base, age_top, timestepsize):
    
    #Nsteps, timestepsize_model = setTimeSteps(age_base, age_top, timestepsize)
    
    # determine time step length
    duration = (age_base - age_top)*1e6*(365*24*60*60)
    Nsteps = np.round(duration/timestepsize)
    
    # minimum timesteps per tectonic time period = 5:
    Nsteps[np.where(Nsteps<5)[0]] = 5
    timestepsize_model = duration/Nsteps
    
    return Nsteps, timestepsize_model


def finiteVolumeFunc(initial_temperatures, diffusivity_s,
    heatProduction_s, cellSize_s, timestepsize,
    Nsteps, heatFlow, surfaceTemperature):
    
    Nstrat = len(initial_temperatures)
    
    # set up finite - volume mesh
    nx = len(initial_temperatures)
    # set up irregular 1D mesh
    mesh = fipy.Grid1D(nx = len(cellSize_s), dx = cellSize_s)       
    # set up temperature variable
    variable = fipy.CellVariable(name = 'T', mesh = mesh,
                                value = initial_temperatures)   
        
    # use harmonic mean of cell based diffusivity for cell - faces:
    diffusivity_cellFace = np.zeros((Nstrat + 1))
    diffusivity_cellFace[1: - 1] = 2. / (1./diffusivity_s[1:] +\
                                        1./diffusivity_s[: - 1])
    diffusivity_cellFace[0] = diffusivity_s[1]
    diffusivity_cellFace[-1] = diffusivity_s[ - 2]
    
    # set diffusivity
    D = fipy.FaceVariable(mesh = mesh, value = diffusivity_cellFace)
    
    # set boundary conditions
    boundaryCondition_top = fipy.FixedValue(
                                        faces = mesh.getFacesLeft(),
                                        value = surfaceTemperature)
    boundaryCondition_bottom = fipy.FixedFlux(
                                        faces = mesh.getFacesRight(),
                                        value = heatFlow)
    
    # add source term:
    SourceTerm = fipy.variables.cellVariable.CellVariable(mesh=mesh,
                                            value = heatProduction_s)
    
    # initialize equation
    diffTerm = fipy.ImplicitDiffusionTerm(coeff = D)
    diffEqTransient = diffTerm+SourceTerm == fipy.terms.TransientTerm()    
    
    # solve finite volume equation for each timestep:
    for step in range(int(Nsteps)):
        diffEqTransient.solve(  var = variable, 
                                boundaryConditions =
                                (boundaryCondition_top,
                                boundaryCondition_bottom),
                                dt = timestepsize)
    
    tempOut = np.asarray(variable[:])
    
    return tempOut



def simulateTemperatures(Nsteps, timestepsize_model, cellSize,
    diffusivity, Kbulk, surfaceTemperature, heatFlow_model,
    heatFlow_original, heatProduction, subsidence_code,
    useSteadyStateHeatFlow):
    
    Ntectonic_steps, Nstrat = cellSize.shape
    temperature = np.zeros((Ntectonic_steps, Nstrat))
    
    temperature[0, 0] = surfaceTemperature[0]
    
    totalTime = 0
    for i in xrange(1, Ntectonic_steps):
        if subsidence_code[i] == 'subsidence' or\
                subsidence_code[i] == 'add_formation':
            temperature[i, 1:] = temperature[i - 1, : - 1]
            temperature[i, 0] = surfaceTemperature[i]
        elif subsidence_code[i] == 'erosion':
            temperature[i, : - 1] = temperature[i - 1, 1:]
    
        
        # find lowest unit
        lowest_unit = np.where(cellSize[i] == 0)[0][0]
        
        if lowest_unit > 2:
            
            # calculate temperature evolution
            if useSteadyStateHeatFlow == False:
                # calculate using implicit finite volume function:
                temperature_output = finiteVolumeFunc(
                    temperature[i, :lowest_unit],
                    diffusivity[i, :lowest_unit],
                    heatProduction[i, :lowest_unit],
                    cellSize[i, :lowest_unit], timestepsize_model[i],
                    Nsteps[i], heatFlow_model[i], surfaceTemperature[i])
            else:
                # calculate using analytical steady - state solution:
                cellDepths = calculateDepths_mid(
                            cellSize[i,:lowest_unit] ) +\
                            cellSize[i, :lowest_unit]*0.5
                temperature_output = calculateSteadyStateTemperatures(
                    heatFlow_original[i], surfaceTemperature[i],
                    cellDepths, Kbulk[i, :lowest_unit])
            
            
            temperature[i, :lowest_unit] = temperature_output.copy()
        
        totalTime += Nsteps[i]*timestepsize_model[i]/(1e6*365*24*60*60)
        
        print ('%s/%s, %i , %0.3f My, T =  %0.1f - %0.1f C, D = %0.2e,ST = %0.1f, HF = %0.2e'
                %(i + 1, Ntectonic_steps, lowest_unit, totalTime,
                    temperature[i, :lowest_unit].min(),
                    temperature[i, :lowest_unit].max(),
                    diffusivity[i, :lowest_unit].mean(),
                    surfaceTemperature[i], heatFlow_model[i]))

    return temperature
    
    
def transformTemperatureArray(temperature, subsidence_code):
    
    Ntectonic_steps, Nstrat = temperature.shape
    sampleTemperatures = np.zeros(temperature.shape)
    NdepthSteps = 0
    column_top = 0
    column_lowest = np.where(temperature[-1, :]>0)[0][-1]
    Nstrat = column_lowest
    for i in xrange(Nstrat + 1):
        loc = i ;erosionSwitch = False
        for j in xrange(Ntectonic_steps - 1, 0, - 1):
            if loc<0:
                erosionSwitch = True
            
            if loc >=  0 and erosionSwitch == False:
                sampleTemperatures[j, i] = temperature[j, loc].copy()
            if subsidence_code[j] == 'subsidence' or subsidence_code[j] == 'add_formation':                
                # insert new stratigraphic unit and surface temperature
                loc  -=  1 
            elif subsidence_code[j] == 'erosion':
                loc += 1
            
    return sampleTemperatures
    

def plotThermalParameters(cellSize, Kbulk, density, diffusivity,
                            age_base):
    
    pl.clf()
    print '\tcreating fig of thermal params'
    cellDepths_mid = calculateDepths_mid(cellSize)
    maxAge = age_base.max()
    maxDepth = cellDepths_mid.max()
    
    for i in xrange(3):
        pl.subplot(2, 2, i + 1)
        if i == 0:
            z = Kbulk
            title_ = 'Thermal conductivity'
        elif i == 1:
            z = density
            title_ = 'Density'
        elif i == 2:
            z = np.log10(diffusivity)
            title_ = 'Log Diffusivity'
        xt = age_base
        x = np.zeros(cellDepths_mid.shape)
        a, b = cellDepths_mid.shape
        for j in xrange(b):
            x[:, j] = xt
        y = cellDepths_mid
        
        zm = np.ma.masked_where(cellSize == 0, z)
        xm = np.ma.masked_where(cellSize == 0, x)
        ym = np.ma.masked_where(cellSize == 0, y)
    
        #zi = griddata(xm, ym, zm, xi, yi)
    
        pl.contourf(xm, ym, zm)
        pl.colorbar()
        pl.ylim(maxDepth, 0)
        pl.title(title_)
    pl.savefig('thermal_params_debug.png')
    print '\tsaved debug figure thermal_params_debug.png'
    #pl.show()
    
    return
    


def simulateBurialAndTemperatureHistory(subsidence_code, 
    stratCodesTable, strat_units, age_base, age_top, 
    cellSize_initial, Phi0, Phic, NcompactionIterations, 
    heatflowMarkers, heatflowHist, surfaceTempHistMarkers, 
    surfaceTempHist, heatProduction, useSteadyStateHeatFlow, 
    changeThickness_units, changeThickness_value, 
    changeThickness_timing, 
    Kclay, Ksand, Kchalk, Korg, Kevap, Kwater, 
    density_water, density_formation, c_water, c_formation, 
    timestepsize, 
    useVR = True):

    print 'arrange stratigraphic units and cell sizes'
    # get stratigraphic units
    if verbose == True:
        print 'test burial sequence'
        checkBurial(strat_units, subsidence_code, 
                    cellSize_initial, cellSize_initial)
    strat_steps = simulateStratUnits(strat_units, subsidence_code)
    
    # get first estimate of cell size
    Nstrat = len(strat_steps)*2
    cellSize = simulateBurial(cellSize_initial, subsidence_code, Nstrat)

    # get initial cell sizes:
    cellSize_initial_steps = simulateBurial(cellSize_initial, 
                                            subsidence_code, Nstrat)

    #simulation_time = calculateSimulationTime(age_base[50:], age_top[50:])
    #print simulation_time
    #print age_base[50:]
    #exit()

    # add fault movement:
    if changeThickness_units !=  []:
        print 'adding fault movement'
        cellSize_initial_steps = changeCellThickness(
            cellSize_initial_steps, strat_units, strat_steps, 
            changeThickness_timing, changeThickness_units, 
            changeThickness_value)

    # find lithologies
    print 'finding lithologies'
    fractionArray = findLithologies(Nstrat, strat_steps,
                                        stratCodesTable)
    
    # calculate compacted cell thickness and porosity
    print 'calculating compaction'
    porosity, cellSize = simulateCompaction(cellSize_initial_steps,
                                            fractionArray, Phi0, Phic)

    # set thermal parameters
    print 'setting thermal parameters'
    Kbulk = calculateThermalConductivity(   fractionArray, porosity, 
                                            Kclay, Ksand, Kchalk, Korg,
                                            Kevap, Kwater )
    density = (porosity*density_water) + (-porosity + 1) *\
                            density_formation  # set density
    heatCapacity = (porosity*c_water) + (-porosity + 1) *\
                            c_formation # set specific heat
    diffusivity = Kbulk/(density*heatCapacity)
    
    if verbose == True:
        plotThermalParameters(cellSize, Kbulk, density,
                                diffusivity, age_base)
    
    # set heat production
    Ntectonic_steps = len(subsidence_code)
    heatProduction_model = np.ones((Ntectonic_steps, Nstrat)) *\
                                                        heatProduction

    # set temperature boundary conditions:
    print 'set boundary conditions'        
    surfaceTemperature = setSurfaceTemp(surfaceTempHist,
                                surfaceTempHistMarkers, strat_units)
    heatFlow = setHeatFlow(heatflowMarkers, heatflowHist, strat_units)

    # correct heat flow and heat production
    # to adhere to FiPy unit definition:
    heatProduction_model_adj =\
                    heatProduction_model / (heatCapacity*density)
    #print heatFlow
    #exit()
    heatFlow_model = np.zeros(heatFlow.shape)
    for i in xrange(Ntectonic_steps):
        lowest_unit = np.where(cellSize[i] == 0)[0][0] - 1
        heatFlow_model[i] = heatFlow[i] / ( heatCapacity[i, lowest_unit]
                                            *density[i, lowest_unit] )
        #heatFlow_model[i] = heatFlow[i]
        #print heatCapacity[i, lowest_unit - 1]
        #print heatCapacity[i, :lowest_unit].mean()
        
        #print density[-1, lowest_unit]
    
    # set time steps
    Nsteps, timestepsize_model = setTimeSteps( age_base, age_top,
                                                timestepsize )
        
    # simulate temperature evolution
    print 'simulating temperature evolution'
    temperature = simulateTemperatures( Nsteps, timestepsize_model,
                                        cellSize, diffusivity, Kbulk,
                                        surfaceTemperature,
                                        heatFlow_model, heatFlow,
                                        heatProduction_model_adj,
                                        subsidence_code,
                                        useSteadyStateHeatFlow)
    
    return  temperature, cellSize, strat_steps, surfaceTemperature,\
            heatFlow, heatProduction, diffusivity, porosity


def calculateResidual_temperature(sampleTemperatures, sampleDepth_,
    DSTgaugeDepth, DSTintervalDepthTop, DSTintervalDepthBottom,
    DSTtemperature, logDepth, logTemp,
    BHTcDepth, BHTcTemp, BHTdepth, BHTtemp, doProb):

    print 'calculate model fit to temperature data'
    
    ###############################
    # calculate difference observed and simulated teperature data
    ###################################
    P_allScenarios = 0 ; Ntotal = 0
    
    #####################################
    # unCorrected bottom hole temperatures
    #####################################
    uncRange = 5. ; maxRange = 20.
    Ndata = len(BHTtemp) ; P_BHT = np.zeros((Ndata))
    #print '%s BHTc data' %Ndata
    
    if Ndata>0:
        for i in xrange(Ndata):
            simulated = linearInterpolation(BHTdepth[i], sampleDepth_,
                                                    sampleTemperatures)
            
            if simulated<BHTtemp[i]:
                difference_ = BHTtemp[i] - simulated
            # simulated temperature >x degrees higher than raw bht
            elif (simulated - BHTtemp[i])>maxRange:     
                difference_ = BHTtemp[i] + maxRange - simulated
            else:
                difference_ = 0
            if verbose == True:
                print 'raw BHT = %0.1f, sim = %0.1f, diff = %0.1f'\
                    %(BHTtemp[i], simulated, difference_)
            #SE = residual**2.*weight
            # normalized difference =  1 when simulated AFT equals 
            #  95%  of confidence interval
            difference_norm = abs(difference_/uncRange)*1.96   
            if doProb == True:
                P_BHT[i] = (1 - stats.norm.cdf(abs(difference_norm)))*2
            else:
                P_BHT[i] = difference_norm**2
            P_allScenarios = P_allScenarios + P_BHT[i]
            Ntotal += 1
        if verbose == True:
            print '%s raw BHT data, P = ' %(Ndata)
            print P_BHT
            
    else:
        P_BHT = np.zeros((1))
    
    #####################################
    # corrected bottom hole temperatures
    #####################################
    uncRange = 5. ; Ndata = len(BHTcTemp)
    P_BHTc = np.zeros((Ndata))
    #print '%s BHTc data' %Ndata
    if Ndata>0:
        for i in xrange(Ndata):
            simulated = linearInterpolation(BHTcDepth[i], sampleDepth_,
                                                    sampleTemperatures)
            difference_ = BHTcTemp[i] - simulated
            #SE = residual**2.*weight
            difference_norm = abs(difference_/uncRange)*1.96
            if doProb == True:
                P_BHTc[i] = (1 - stats.norm.cdf(abs(difference_norm)))*2
            else:
                P_BHTc[i] = difference_norm**2

            P_allScenarios = P_allScenarios + P_BHTc[i]
            Ntotal += 1
        if verbose == True: print '%s corrected BHT data, P = ' %(Ndata)
    else:
        P_BHTc = np.zeros((1))
         
    #####################
    # DST's
    ######################
    uncRange = 3.
    
    Ndata = len(DSTintervalDepthTop) 
    ybest = (DSTintervalDepthTop + DSTintervalDepthBottom)/2
    x = DSTtemperature
    P_DST = np.zeros((Ndata))
    #print '%s DST data' %Ndata        
    if Ndata>0:
        for i in xrange(Ndata):
            simulated = linearInterpolation(ybest[i],
                                    sampleDepth_, sampleTemperatures)
            difference_ = x[i] - simulated
            # i.e. normalized difference =  1 when simulated AFT equals 
            #  95% confidence interval
            difference_norm = abs(difference_/uncRange)*1.96   
            if doProb == True:
                P_DST[i] = (1 - stats.norm.cdf(abs(difference_norm)))*2
            else:
                P_DST[i] = difference_norm**2
            
            P_allScenarios = P_allScenarios + P_DST[i]
            Ntotal += 1
        if verbose == True: print '%s DST data, P = ' %(Ndata)
                
    else:
        P_DST = np.zeros((1))

    if Ntotal !=  0:
        if doProb == False:
            residual = math.sqrt((P_allScenarios/Ntotal))
        else:
            residual = (P_allScenarios/Ntotal)
    else:
        residual = 0

    if verbose == True:
        print 'error DSTs:', P_DST
        print residual
        print doProb
        print P_BHTc
    
    return residual


def calculateResidual_VR(VRDepth, VRValue,
    VRValue_std, cellDepths_present,
    VR_simulated, doProb):
    # calculate difference between observed (VRValue) and simulated (VR_simulated) VR values 
    
    SE = 0 ; N = 0
    P = np.zeros((len(VRDepth)))
    for i in xrange(len(VRDepth)):
        if VRValue_std[i]>0 and VRValue_std[i]<100.:
            uncRange = VRValue_std[i]*1.96
        else:
            uncRange = 0.1
        VR_sim_int = linearInterpolation(VRDepth[i], cellDepths_present,
                                            VR_simulated)
        difference_ = VRValue[i] - VR_sim_int
        difference_norm = abs(difference_/uncRange)*1.96   # i.e. normalized difference =  1 when simulated AFT equals 95% confidence interval
        if doProb == True:
            P[i] = (1 - stats.norm.cdf(abs(difference_norm)))*2
        else:
            P[i] = difference_norm**2
        
    #residual = 1 - P.mean()
    residual = math.sqrt(P.mean())
    
    return residual


def addProvenanceAges(provenance_age, provenance_age_step,
                                    AFTtime_, AFTtemp_):
    
    # same as first but always uses a linear change in temperature ie. temperature at 2nd timestep is 60 degrees, time can be calibrated
    
    ################################################################################
    # add provenance (pre-deposition) temperature history to time + temperature arrays
    ################################################################################

    AFTtime_prov_ = np.zeros((len(AFTtime_) + 2))
    AFTtemp_prov_ = np.zeros((len(AFTtime_) + 2))
    
    AFTtime_prov_[2:] = AFTtime_ ; AFTtemp_prov_[2:] = AFTtemp_

    if provenance_age_step == 0: doLinear = True
    else: doLinear = False
    
    # check if provenance age is higher than stratigraphic age
    if provenance_age_step<AFTtime_[0]:
        provenance_age_step = (AFTtime_[0]) + 1.        
    if provenance_age<AFTtime_[0] or provenance_age<provenance_age_step:
        provenance_age = provenance_age_step + 1.
        
    ######################################
    # add provenance age to time array
    ######################################
    AFTtime_prov_[0] = provenance_age 
    AFTtime_prov_[1] = provenance_age_step
    
    #######################################################
    # set temperature during provenance temperature history
    #######################################################
    AFTtemp_prov_[0] = 120. # starting temp at provenance age 
    AFTtemp_prov_[1] = AFTtemp_prov_[2]    # surface temperature at 2nd provenance time step
    
    return AFTtime_prov_, AFTtemp_prov_


def wrapProvenanceParams(provenance_ages,
                                    provenance_age_step, NAFTsamples):
    #print NAFTsamples
    #global Ncomponents
    dummy, Ncomponents = np.shape(provenance_ages)
    parameters = np.zeros((Ncomponents*2*NAFTsamples))
    parameters[:NAFTsamples*Ncomponents] = provenance_ages[:NAFTsamples].flatten() 
    parameters[NAFTsamples*Ncomponents:] = provenance_age_step[:NAFTsamples].flatten()
    
    return parameters


def unwrapProvenanceParams(parameters, NAFTsamples):
    
    Ncomponents = len(parameters)/NAFTsamples/2.
    provenance_ages = np.zeros((NAFTsamples, Ncomponents))
    provenance_age_step = np.zeros((NAFTsamples, Ncomponents))
    
    Nprov = Ncomponents*NAFTsamples
    for sampleNo in xrange(NAFTsamples):
        provenance_ages[sampleNo, :] = parameters[sampleNo*Ncomponents:(sampleNo + 1)*Ncomponents]
        provenance_age_step[sampleNo, :] = parameters[Nprov + sampleNo*Ncomponents:Nprov + (sampleNo + 1)*Ncomponents]
        
    return provenance_ages, provenance_age_step


def evaluateProvenanceHistoryFit(parameters, useCalibrate, NAFTsamples, 
    AFTtime_list, AFTtemp_list, kineticParam, kineticValues,
    binsize, useCaxis, AFTages, AFTages_min, AFTages_max,
    trackLengths, useModelFit):

    Ncomponents = len(parameters)/NAFTsamples/2
           
    Nprovenance_ages = Ncomponents
    
    a,N_compositions = kineticValues.shape
    
    provenance_ages, provenance_age_step =\
                    unwrapProvenanceParams(parameters, NAFTsamples)
    
    ###############################################################
    # simulate apatite fission track ages and length distributions
    ###############################################################
    AFTtime_prov_all = []
    AFTtemp_prov_all = []
    
    trackLengthPDF_ = np.zeros((NAFTsamples, int(20/binsize), Nprovenance_ages, N_compositions))
    simulatedAFTages_ = np.zeros((NAFTsamples, Nprovenance_ages, N_compositions))
    age_reduction_prov = np.zeros((NAFTsamples, Nprovenance_ages, N_compositions))
    age_reduction_basin = np.zeros((NAFTsamples, Nprovenance_ages, N_compositions))
    
    if verbose == True: fig = pl.figure()

    #count = 0
    for sampleNo in xrange(NAFTsamples):
        
        AFTtime_ = AFTtime_list[sampleNo]   
        AFTtemp_ = AFTtemp_list[sampleNo]

        AFTtime_prov = np.zeros((len(AFTtime_) + 2, Nprovenance_ages))
        AFTtemp_prov = np.zeros((len(AFTtemp_) + 2, Nprovenance_ages))
                
        #############################################################
        # calculate apatite fission track length distribution and age
        #############################################################  
        count = 0  
        for j in xrange(Nprovenance_ages):
                
            l = j
            
            # add provenance temperature history:
            AFTtime_prov_, AFTtemp_prov[:, j]\
             = addProvenanceAges(provenance_ages[sampleNo, j],
                    provenance_age_step[sampleNo, j],
                    AFTtime_, AFTtemp_)
            #  change from yr BP to My after start of 
            #  temperature history
            AFTtime_prov[:, j] = -AFTtime_prov_ + AFTtime_prov_.max() 
            
            # add compositions:
            for k in xrange(N_compositions):        
                
                print '-'*10
                print 'sample %s, prov scen %i, %s=%0.2f'\
                        %((sampleNo+1), (j+1), kineticParam[sampleNo],(kineticValues[sampleNo,k]))
                # calculate AFT age and track length distribution
                trackLengthPDF_[sampleNo, :, j, k],\
                simulatedAFTages_[sampleNo, j, k],\
                l_mean, l_mean_std, rm, rc, r_cmod, rho_age, dt\
                 = AFTannealingLib.simulate_AFT_annealing( 
                            AFTtime_prov[:, j], AFTtemp_prov[:, j],
                            kineticValues[sampleNo,k],
                            kineticParameter=kineticParam[sampleNo])
                
                    # old function
                    #trackLengthPDF_[sampleNo, :, j, k],\
                    #simulatedAFTages_[sampleNo, j, k],\
                    #l_mean, l_mean_std, rm, rc, r_cmod, rho_age, dt\
                    #= AFTlib_old.calculate_AFT_Ketcham(
                    #    AFTtime_prov[:, j], AFTtemp_prov[:, j],
                    #    'Ketcham2007', -999, kineticValues[k],
                    #    -999, 0.25, False )
                
#                except:
#                    print 'error calculating AFT data'
#                    exit()
                #pdb.set_trace()
        # store simulated temperature history in list:
        AFTtime_prov_all.append(AFTtime_prov)
        AFTtemp_prov_all.append(AFTtemp_prov)
    
    #exit()
    
    if useModelFit == True:
        ###########################
        # calculate residuals
        ###########################
        doProb = True
        
        ##########################################################
        # calculate residual of simulated vs measured AFT ages:
        ##########################################################
        if useCalibrate == False:
            GOF_AFTages = calculateResidual_AFTages( AFTages, 
                            AFTages_min, AFTages_max,
                            simulatedAFTages_, doProb)
        RMSE_AFTages = calculateResidual_AFTages(AFTages,
                            AFTages_min, AFTages_max,
                            simulatedAFTages_, False)
        
        #AFTages, AFTages_min, AFTages_max, trackLengths, binsize
        
        
        ##########################################################
        # calculate residual of simulated vs measured AFT lengths:
        ##########################################################
        if useCalibrate == False:
            GOF_AFTlengths = calculateResidual_AFTlengths(
                                    trackLengthPDF_, trackLengths,
                                    doProb, binsize )
                                    
        RMSE_AFTlengths, TLSD_difference, MTL_observed,\
        MTL_simulated_min, MTL_simulated_max\
         = calculateDifferenceMeanTrackLengths(trackLengths,
                                            trackLengthPDF_, binsize)
        
        print '\tRMSE FT lengths: %0.2f' %RMSE_AFTlengths
        
        # total residual:
        residual = (RMSE_AFTages + (RMSE_AFTlengths*10))
    
        if useCalibrate == True:
            print '\n' + ' + '*int(residual*4)
            print 'residual =  %0.3f, components %s' %(residual, Ncomponents)
            print ' + '*10
            return residual
        else:
            return AFTtime_prov_all, AFTtemp_prov_all,\
            trackLengthPDF_, simulatedAFTages_, GOF_AFTages,\
            GOF_AFTlengths, RMSE_AFTages, RMSE_AFTlengths
    else:
        # no model fit, return only simulated fission track data:
        return AFTtime_prov_all, AFTtemp_prov_all,\
                trackLengthPDF_, simulatedAFTages_
                

def showModelCellSizes(strat_steps, wellStratCode, cellSize_present,
    cellSize, wellStratDepth_baseTVD, wellStratDepth_topTVD):
    
    # show comparison of model and well strat. units thickness
    th = 0 ; th_present_model = 0 ; count = 0
    for i in xrange(0, len(strat_steps[-1])):
        s_mod = strat_steps[-1][i]
        if '-' in s_mod: s_mod = s_mod[:s_mod.index('-')]
        s_well = wellStratCode[count]
        if s_mod !=  s_well:
            th_well = (wellStratDepth_baseTVD[count] - wellStratDepth_topTVD[count])
            print '%s - %s/%s thickness, model =  %0.0f, well = %0.0f / %0.0f'\
             %(count, s_mod_old, s_well, th, th_well, th_present_model)
            count += 1 ; th = cellSize[-1, i] ; th_present_model = cellSize_present[i]
        else:
            th = th + cellSize[-1, i] ; th_present_model += cellSize_present[i]
        s_mod_old = s_mod
    
    return


def showCellSizes_present( cellSize_present_all, cellSize_all,
                            strat_steps ):
    
    Nstrat = len(strat_steps)
    for i in xrange(Nstrat):
        print '%s, cellsize model: %0.1f, cellsize well: %0.1f '\
        %(strat_steps[i], cellSize_all[i], cellSize_present_all[i])
    print len(strat_steps)
    print len(cellSize_all)
    print len(cellSize_present_all)
    
    return 


def runModelScenario(
    maxDuration, maxCellSize,
    lateCret_erosion_start, lateCret_erosion_end,
    lateCret_erosion,
    prePerm_erosion_start, prePerm_erosion_end, prePermEroded,
    erosion_stratCode, erosion_stratCode2,
    erosion_start, erosion_end, eroded_thicknes, 
    heatflowMarkers, heatflowHist,
    surfaceTempHistMarkers, surfaceTempHist, heatProduction,
    useSteadyStateHeatFlow, calibrateProvenanceScenarios,
    provenance_ages, provenance_age_step, BHTDepth, BHTValue,
    DSTgaugeDepth, DSTintervalDepthTop, DSTintervalDepthBottom,
    DSTtemperature, logDepth, logTemp, BHTcDepth, BHTcTemp,
    VRDepth, VRValue, VRValue_std,
    AFTsample_id, AFTDepth_base, AFTDepth_top,
    trackLengths, trackAngles,
    AFTages, AFTages_min, AFTages_max, binsize, 
    h5file, Phi0, Phic, VR_timestep,
    changeThickness_units, changeThickness_value,
    changeThickness_timing, maxATthickness_, 
    well , eroded_thickness,
    kineticParam, kineticValues,
    Kclay, Ksand, Kchalk, Korg, Kevap, Kwater,
    density_water, density_formation, c_water, c_formation,
    timestepsize,
    useVR = True, useCaxis = False ):

    
    NAFTsamples = len(AFTsample_id)
    dummy, Nprovenance_ages = np.shape(provenance_ages)
    a, N_compositions = kineticValues.shape
    
    trackLengthPDF_ = np.zeros((NAFTsamples, int(20/binsize),
                                Nprovenance_ages, N_compositions))
    simulatedAFTages_ = np.zeros((NAFTsamples, Nprovenance_ages,
                                N_compositions))
    
    print 'start preparing well data'
    cellSize_initial, cellSize_present, strat_units,\
    subsidence_code, age_base, age_top,\
    stratCodesTable, wellStratCode, wellStratDepth_topTVD,\
    wellStratDepth_baseTVD, prePerm,\
    totalATerosion, totalSLerosion\
     = prepareWellData(
        h5file, well, maxDuration, maxCellSize,
        lateCret_erosion_start, lateCret_erosion_end, lateCret_erosion,
        erosion_stratCode, erosion_stratCode2,
        erosion_start, erosion_end, eroded_thickness,
        prePerm_erosion_start, prePerm_erosion_end, prePermEroded,
        Phi0, Phic, maxATthickness_ )
    print 'well data read from hdf5 database'

    # check duration:
    duration = age_base - age_top
    if duration.min() <=  0:
        print 'error preparing well subsidence data, <0 duration timestep'
        print duration
        sys.exit()

    
    # remove all DST data deeper than deepest strat unit
    select_ = np.where(DSTintervalDepthBottom<wellStratDepth_baseTVD[-1])[0]
    DSTgaugeDepth = DSTgaugeDepth[select_]
    DSTintervalDepthTop = DSTintervalDepthTop[select_]
    DSTintervalDepthBottom = DSTintervalDepthBottom[select_]
    DSTtemperature = DSTtemperature[select_]
    
    ######################################################################################
    # remove all episodes older than oldest formation in well from heat flow history array
    ######################################################################################
    delete_index = 999
    for i in xrange(len(heatflowHist)):
        HFstrat = heatflowMarkers[i] ; a = len(HFstrat)
        doDelete = True
        for j in xrange(len(strat_units)):
            if strat_units[j][:a] == HFstrat or heatflowMarkers[i] == 'all': doDelete = False
        if doDelete == True:
            delete_index = i

    if delete_index !=  999:
        print 'deleted marker %s from heat flow history, not present in well strat' %heatflowMarkers[delete_index]
        heatflowMarkers = np.delete(heatflowMarkers, delete_index)
        heatflowHist = np.delete(heatflowHist, delete_index)
    
    # make figure of subsidence
    if verbose == True:
        print 'creating figure test_ages.png of subsidence:'
        pl.figure()
        pl.subplot(2, 1, 1)
        pl.plot(age_base)
        pl.plot(age_top)
        duration_test = age_base - age_top
        pl.subplot(2, 1, 2)
        pl.plot(duration_test)
        pl.savefig('test_ages.png')
        pl.clf()
            
    ###########################################################
    # finite volume heat flow calculations:
    ###########################################################
    NcompactionIterations = 5
    
    temperature, cellSize, strat_steps, surfaceTemperature, \
    heatFlow, heatProduction_model, diffusivity, porosity\
    = simulateBurialAndTemperatureHistory(subsidence_code
        , stratCodesTable, strat_units, age_base, age_top
        , cellSize_initial, Phi0, Phic, NcompactionIterations
        , heatflowMarkers, heatflowHist, surfaceTempHistMarkers
        , surfaceTempHist, heatProduction, useSteadyStateHeatFlow
        , changeThickness_units, changeThickness_value
        , changeThickness_timing
        , Kclay, Ksand, Kchalk, Korg, Kevap, Kwater
        , density_water, density_formation, c_water, c_formation
        , timestepsize)
    # convert cell sizes found in well to same format as model output
    Nstrat = len(strat_steps)*2
    cellSize_present_steps = simulateBurial(cellSize_present, subsidence_code, Nstrat)

    cellDepths = calculateDepths_mid(cellSize)
    cellDepths_mid = calculateDepths_mid(cellSize)
    cellDepths_base = calculateDepths_base(cellSize)
    cellDepths_present_well = wellStratDepth_baseTVD.copy()
    cellDepths_present_mid = calculateDepths_mid(cellSize_present_steps)
    print 'calculated cell depths: %0.0f' %(cellDepths_base[-1, - 1])
    print 'cell depths from present day well strat horizon depths: %0.0f' %cellDepths_present_well[-1]
    
    if verbose == True or a == 1:
        showModelCellSizes(strat_steps, wellStratCode, cellSize_present, cellSize, wellStratDepth_baseTVD, wellStratDepth_topTVD)
            
    #########################################################
    # convert 2d temperature profiles array to 
    # 2d temperature history array of present day formations
    #########################################################
    
    sampleTemperatures = transformTemperatureArray(temperature, subsidence_code)
    
    if verbose:
        print 'saving file with temperature histories: sampleTemp_debug.csv'
        np.savetxt('sampleTemp_debug.csv', sampleTemperatures)
        

    ####################################################
    # calculate vritinite reflectance:
    ####################################################
    
    if useVR == True:
        # calculate simulated VR
        # determine number of formations/grid_cells at present
        Nsamples = len(np.nonzero(sampleTemperatures[-1, :])[0])
        #global age_base
        VR_simulated = calculateVR(sampleTemperatures, subsidence_code, Nsamples, age_base, VR_timestep)
            
    #####################################################################
    # extract temperature history apatite samples from simulation results
    #####################################################################
    print 'calculating temperature history fission track samples'
    AFTtime_list = [] ; AFTtemp_list = []
    NAFTsamples = len(AFTDepth_base)
    for sampleNo in xrange(NAFTsamples):
        # calculate simulated time - temperature path
        AFTtime, AFTtemp = calculateSimulatedAFT_precise(AFTDepth_base[sampleNo], AFTDepth_top[sampleNo]\
        , sampleTemperatures, age_base, age_top, cellDepths_present_mid[-1], surfaceTemperature, strat_steps, subsidence_code)
        AFTtime_list.append(AFTtime)
        AFTtemp_list.append(AFTtemp)

    
    
    
    ###############################################################
    # simulate apatite fission track ages and length distributions
    ###############################################################
    test = False
    
    dummy, Ncomponents = np.shape(provenance_ages)
    min_residual = 0.7 ; residual = 0 ; Nsteps = 0
    
    # set calibration parameters
    parameters = wrapProvenanceParams(provenance_ages, 
                                    provenance_age_step, NAFTsamples)
    if calibrateProvenanceScenarios == True:
        print 'calibrating provenance thermal history:'
        print '\tinitial parameters: ' 
        print '\t', parameters
    
    ##########################################
    # calibrate provenance temperature history
    ##########################################
    if calibrateProvenanceScenarios == True: 
        opt_results\
         = opt.fmin(evaluateProvenanceHistoryFit, parameters,
            args = (True, NAFTsamples, AFTtime_list, AFTtemp_list,
            kineticParam, kineticValues, binsize, useCaxis,
            AFTages, AFTages_min, AFTages_max,
            trackLengths, True), ftol = 0.01, xtol = 0.01,
            disp = True, full_output = True, retall = True)
        
        # read final parameters
        parameters = opt_results[0]
        print '\tfinal calibrated parameters: ', parameters
        
        # read provenance age parameteres from calibrated parameters
        provenance_ages, provenance_age_step  =\
            unwrapProvenanceParams(parameters, NAFTsamples)
        
        print '\tprovenance ages:'
        print '\t', provenance_ages
        
        print 'provenance_age_step:'
        print '\t', provenance_age_step
            
    # check if residual is low enough:
    AFTtime_prov_all, AFTtemp_prov_all, trackLengthPDF_,\
    simulatedAFTages_, GOF_AFTages, GOF_AFTlengths,\
    RMSE_AFTages, RMSE_AFTlengths\
     = evaluateProvenanceHistoryFit(parameters, False, NAFTsamples,
            AFTtime_list, AFTtemp_list, kineticParam, kineticValues,
            binsize, useCaxis, AFTages, AFTages_min, AFTages_max,
            trackLengths, True)
    
    residual = GOF_AFTages
    
    ###########################
    # calculate residuals
    ###########################
    doProb = True
    GOF_temperature = calculateResidual_temperature(sampleTemperatures[-1, :Nsamples].copy(), cellDepths[-1, :Nsamples]\
    , DSTgaugeDepth, DSTintervalDepthTop, DSTintervalDepthBottom, \
    DSTtemperature, logDepth, logTemp, BHTcDepth, BHTcTemp, BHTDepth, BHTValue, doProb)
    RMSE_temp = calculateResidual_temperature(sampleTemperatures[-1, :Nsamples].copy(), cellDepths[-1, :Nsamples]\
    , DSTgaugeDepth, DSTintervalDepthTop, DSTintervalDepthBottom, \
    DSTtemperature, logDepth, logTemp, BHTcDepth, BHTcTemp, BHTDepth, BHTValue, False)
              
    ##########################################################
    # calculate VR residual (observed - simulated values)
    ##########################################################
    GOF_VR = calculateResidual_VR(VRDepth, VRValue, VRValue_std\
    , cellDepths[-1, :Nsamples], VR_simulated, doProb)
    RMSE_VR = calculateResidual_VR(VRDepth, VRValue, VRValue_std\
    , cellDepths[-1, :Nsamples], VR_simulated, False)
    
    print '-'*10
    print 'goodness of fit:'
    print 'temperature: %0.2f' %GOF_temperature
    print 'VR: %0.2f' %GOF_VR
    print 'AFT ages avg: %0.2f' %GOF_AFTages.mean()
    for j in xrange(len(GOF_AFTages)):
        print '\tsample %s, age: %0.2f, length: %0.3f'\
                %(j, GOF_AFTages[j], GOF_AFTlengths[j])
    
    
    return temperature, AFTtime_prov_all, AFTtemp_prov_all,\
    AFTtime_list, AFTtemp_list,\
    cellDepths, cellSize, cellSize_present_steps,\
    subsidence_code, strat_units, strat_steps, age_base, age_top,\
    VR_simulated, heatFlow, sampleTemperatures, \
    lateCret_erosion, totalSLerosion, totalATerosion, \
    GOF_temperature, GOF_VR, GOF_AFTlengths, GOF_AFTages, \
    RMSE_temp, RMSE_VR, RMSE_AFTlengths, RMSE_AFTages, \
    AFTages, AFTages_min, AFTages_max, trackLengths,\
    trackAngles, trackLengthPDF_, simulatedAFTages_,\
    provenance_ages, provenance_age_step,\
    diffusivity, porosity
    
def timeDuration(seconds):
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    years, days = divmod(days, 365.242199)
 
    return years, days, hours, minutes, seconds


def constructHeFTyFile(time,temperature,fileName):
    """
    create HeFTy time-temperature input file
    
    input:
    time:           time in Ma bp
    temperature:    temperature in degr. C
    fileName:       name of HeFTy file to write results to
    """    
        
    outputStr = ''
    for i in xrange(len(time)):
        outputStr += '%0.3f, %0.2f' %(time[i],temperature[i])
        if i<(len(time)-1):
            outputStr += '\n'
    
    fout = open(fileName,'w')
    fout.write(outputStr)
    fout.close()
    
    return
