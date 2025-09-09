'''


'''

import itertools, pdb, math
import numpy as np
from numba import jit


#@jit(nopython=True)
def easyRo(times, temperatures_in, vr_method='easyRo', debug=False):
    
    """
    calculate vitrinite reflectance using the easy%Ro algorithm of 
    Sweeney & Burnham (1990) AAPG Bulletin 74(10)
    
    this function has been verified by reproducing the %Ro values 
    provided in Figure 8 of Sweeney & Burnham (1990)
    
    input:
        times           time in My, first item assumed to be 0 My
        temperatures    temperature at each timestep, in degrees C
        vr_method       VR method, choose either 'easyRo' or 'basinRo'. default is 'easyRo'
        
    returns:
        Ro              calculated %Ro values at each timestep
    """

    Myr = 1e6 * 365.25 * 24 * 60 * 60

    timesteps = len(times)
    
    # convert T to Kelvin
    temperatures = np.copy(temperatures_in)+273.15

    # fixed parameters:
    if vr_method is 'easyRo':
        weights = np.array([0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.06,
                            0.04, 0.04, 0.07, 0.06, 0.06, 0.06, 0.05,
                            0.05, 0.04, 0.03, 0.02, 0.02, 0.01 ])
    elif vr_method is 'basinRo':
        weights = np.array([0.0185,	0.0143,	0.0569,	0.0478,	0.0497,	0.0344,	0.0344,
                            0.0322,	0.0282,	0.0062,	0.1155,	0.1041,	0.1023,	0.076,
                            0.0593,	0.0512,	0.0477,	0.0086,	0.0246,	0.0096])

    # activation energy (kcal/mol)
    activationEnergy = np.array([34., 36., 38., 40., 42., 44., 46., 48., 50.,
                                 52., 54., 56., 58., 60., 62., 64., 66., 68.,
                                 70., 72. ])
    
    components = len(activationEnergy)
    
    # preexponential factor (s^-1)
    if vr_method is 'easyRo':
        preExp = 1.0e13
    elif vr_method is 'basinRo':
        preExp = np.exp(60.9856) / Myr
    
    # universal gas constant (cal K-1 mol-1)  
    R = 1.987  
    
    # rearrange time & activation energy arrays
    T_all_components = np.resize(temperatures,(components, timesteps))
    activationEnergies_all_timesteps = np.resize(activationEnergy,(timesteps, components)).T
    
    # ERT
    EdivRT = activationEnergies_all_timesteps*1000.0 / (R*temperatures)
    
    
    # I
    I = preExp * T_all_components * np.exp(-EdivRT) * (1-(EdivRT**2+2.334733*EdivRT+0.250621) / 
                                        (EdivRT**2+3.330657*EdivRT+1.681534) )
    # calculate heating rates

    heatingRates = np.diff(temperatures) / np.diff(times * Myr)
    
    # zero heating rate at first timestep
    heatingRates = np.insert(heatingRates, [0], 0)
    
    # remove zero heating rates
    heatingRates[heatingRates==0] = 1.0e-30
    
    # delta I
    deltaI = np.zeros(I.shape)
    for j in range(1, timesteps, 1):
        deltaI[:, j] = deltaI[:, j-1] + (I[:, j]-I[:, j-1]) / heatingRates[j]
    
    #cumulativeReacted
    cumulativeReacted = weights * (1.0 - np.exp(-deltaI.T))
    for k in range(components):
        cumulativeReacted[deltaI[k, :] > 220.0, k] = weights[k]
    cumulativeReacted[deltaI.T < 1.0e-20] = 0
    
    #sumReacted
    sumReacted = cumulativeReacted.sum(axis=1)
    
    # calculate Ro at each timestep

    if vr_method is 'easyRo':
        Ro = np.exp(-1.6 + 3.7 * sumReacted)
    elif vr_method is 'basinRo':
        Rzero = 0.2104
        Ro = Rzero * np.exp(3.7 * sumReacted)
    
    if debug == True:
        return Ro, sumReacted, cumulativeReacted, deltaI, I, EdivRT

    else:
        return Ro


#@jit(nopython=True)
def easyRo_slow(timeArray, tempArray, debug = False):
    
    """
    calculate vitrinite reflectance using the easy%Ro algorithm of 
    Sweeney&Burnham (1990) AAPG Bulletin 74(10)
    this function has been verified by reproducing the %Ro values 
    provided in Figure 8 of Burnham&Sweeney(1990)
    
    input:    - a 1D array with the time in My (timeArray), 
                    first item assumed to be 0 My
               - a 1D array that contains the temperature at each
                    timestep,  in degrees C (tempArray)
    output: - a 1D array containing the calculated %Ro values at each
                    timestep
                    
    slower, not numpyfied version
    
    way to build this into Rift2D:
    each timestep:
    - model input = temperature
    - calculate heating rate for each grid cell (mesh nodes * 1 memory/operations)
    - calculate I for each component (mesh nodes  * 20 )
    - calculate deltaI for each component (mesh nodes * 20 )
    last timestep:
    - take last deltaI value and calculate cumulative reacted, sum and Ro
    
    """
    
    # fixed parameters:
    weights = np.array([0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.06, 0.04, 0.04\
    , 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.04, 0.03, 0.02, 0.02, 0.01])
    activationEnergy = np.array([34, 36, 38, 40, 42, 44, 46, 48, 50, 52\
    , 54, 56, 58, 60, 62, 64, 66, 68, 70, 72]) # activation energy (kcal/mol)
    components = len(activationEnergy)
    preExp = 1.0e13  # preexponential factor (s^-1)
    R = 1.987  # universal gas constant (cal K-1 mol-1)
    
    # set up arrays to store calculation steps:
    datapoints = len(timeArray)
    EdivRT = np.zeros((components))
    I = np.zeros((components))
    deltaI = np.zeros((components))
    cumulativeReacted = np.zeros((components))
    sumReacted = np.zeros((datapoints))
    Ro = np.zeros((datapoints))
    
    # perform calculation of vitrinite reflectance
    for j in range(datapoints):
        time = timeArray[j]
        T = tempArray[j]
        T = T+273  # convert T to Kelvin
        EdivRT[:] = 0
        I[:] = 0
        deltaI[:] = 0
        cumulativeReacted[:] = 0
        for k in range(components):
            EdivRT[k] = activationEnergy[k]*1000/(R*T)
            ERT = EdivRT[k]
            I_old = I.copy()
            I[k] = preExp*T*math.exp(-ERT)*(1-(ERT**2+2.334733*ERT+0.250621)/(ERT**2+3.330657*ERT+1.681534))
            if j == 0:
                deltaI[k] = 0
            else:
                if j == 0:
                    heatingRate = 1e-30
                else:
                    heatingRate = ((tempArray[j]-tempArray[j-1])/(timeArray[j]-timeArray[j-1]))/(1e6*365*24*60*60)
                    deltaI_old = deltaI.copy()
                if heatingRate == 0:
                    heatingRate = 1e-30
                
                deltaI[k] = deltaI_old[k]+(I[k]-I_old[k])/heatingRate

            if deltaI[k]<1e-20:
                cumulativeReacted[k] = 0
            elif deltaI[k]>220:
                cumulativeReacted[k] = weights[k]
            else:
                cumulativeReacted[k] = weights[k]*(1-math.exp(-deltaI[k]))
        
        sumReacted[j] = np.sum(cumulativeReacted[:])
        Ro[j] = math.exp(-1.6+3.7*sumReacted[j])
        
        #print 'x'*10
        #print j
        #print deltaI
        #print I_old
        #print I
        #print heatingRate
        
        #pdb.set_trace()
        
        if np.isnan(Ro[j]) == True:
            print('!!error in vitrinite function,  no data value at time slice %s: %s' %(j, Ro[j]))
            print('time  =  %s,  temp =  %s' %(time, T))
            print('input time + temperature arrays:')
            print(timeArray)
            print(tempArray)
            print(bla)
    
    if debug == True:
        return Ro, sumReacted, cumulativeReacted, deltaI, I, EdivRT

    else:
        return Ro
