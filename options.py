import numpy as np
import os
c = 2.9979e10

def inputConstants():


    epsilon = 0. #Radiative efficiency
    epsiloneLog = np.log10(0.2)#-0.33#np.log10(2.5e-1) #Fraction of energy in the electrons
    epsilonpLog = -1 #Above, for protons
    epsilone3Log = -3
    epsilonp3Log = -2
    E0 = np.log10(1e53)    #Initial isotropic energy
    nLog = -0 #np.log10(.) #Circumstellar density
    s = 2.0 #CBM power law slope
    Gamma0log = np.log10(150.) #Initial gamma
    eB = -3#np.log10(3.5e-4) #Magnetic microphysical constant in Forward Shock (region 2)
    eB3 = -4   #Magnetic microphysical constant in Reverse Shock (region 3)
    p = 2.1  #Synchrotron emission power law
    logt0 = 1 #Time between R=0 and prompt time   (logarithmic)
    theta0 = np.pi / 180. * 5       #Opening angle
    alpha =  np.pi / 180. * 1          #Off-axis angle
    tN = 1.
    logrphoto = 13.    #10th logarithm of the photosperic radius of the thermal component (cm)
    #E0 = np.log10(8.1e52 / (1-np.cos(theta0)))
    Tprim0 = 5.e5 * 1.60217657e-12   # E [eV]   change first factor only
    N0 = 1.e6
    z = 1.         #Redshift
    WV = 0.    #WV = 0 : open universe    ;   WV = 0.73 : flat universe
    tpromptObs = 10.         #Duration time of prompt phase. Needed for the reverse shock. Input time is observer frame, the module returns the progenitor frame time. Note from 29/1 -14
    promptGamma = 10.       #The lorentz factor during the prompt emission phase. Will be used to transform tpromptObs to tprompt (prog. frame)
    pRS = 2.05          #Ratio between fixed reverse shock and forward shock microphysical ceofficients. Assumes coeffs of forward shock to set coeffs of reverse shock. If you want to use this, set fixedRSFSratio = TRUE in beginning of mainFit.py.
    R0_Clog = 12  #Initial radius of cocoon
    N_Clog = 50      #Number density of cocoon (logarithmic)
    t_outflow_log = 2    #Duration of outflow producing the thermal component (cocoon). Unit seconds in progenitor frame
    Ctheta = 180 /180 * np.pi    #Initial opening angle of cocoon
    GammaClog = np.log10(1.2)             #Logarithm of initial Lorentz factor of cocoon
#freq = [6.5e13,8.2e14]  #Frequency of lightcurve
    return epsilon,epsiloneLog,epsilonpLog,epsilone3Log,epsilonp3Log,E0,nLog,s,Gamma0log,eB,eB3,p,logt0,theta0,alpha,tN,logrphoto,Tprim0,N0,np.log10(tpromptObs/(1+z) / (1-np.sqrt(1-1/promptGamma**2))+ 10**logt0),pRS,R0_Clog,N_Clog,t_outflow_log,Ctheta,GammaClog,z,WV
#np.savetxt('constants.txt',[epsilon,epsiloneLog,epsilonpLog,epsilone3Log,epsilonp3Log,E0,nLog,s,Gamma0,eB,eB3,p,t0,theta0,alpha,tN,logr0,Tprim0,N0,tpromptObs/(1+z) / (1-np.sqrt(1-1/promptGamma**2))+ t0,RSFSratio,z])
#runMainInput = raw_input('Run main programme? ([y]/n): ')
#if (runMainInput == 'y') | (runMainInput == 'Y') | (runMainInput == ''): runMain = True
#else: runMain = False
#if runMain: os.system('ipython mainFit.py')

def userOptions():
    
########################################
#             User options             #
########################################


#Input data
    inputData = '%s/Data_files/'%(os.path.expanduser('~'))  #Where the input data-files are stored. Store the files in a folder under said folder e.g. named 990510fit for GRB 990510. Then specify the burst date to GRBlabel
    fileExtention = 'dat'                                         #File extentions for input data files
    breakName = 'lc.dat'                                            #Type what comes after the frequency in the file name
    inputTimeFormat = 'sec'                                     #What format is input time in? 'day', 'hour', 'min' or 'sec'
    inputFreqUnit = 'Hz'                                        #In what unit the frequency is stated in the filename. 'cm','Hz','GHz','MHz'
    plotInput = True                                              #Print the input data plot to a file? True: yes
    GRBlabel = '990510'                                           #Name of the GRB. If you want manual input, set GRBinput to True
    GRBinput = False

#Options
    allowDataRemoval = False       #True: allow the programme to ask if the user wants to remove previously saved data
    gridReducer = False            #Use grid reducer?
    surfacePlot = True
    
    
    thermalComp = False             #Include a thermal component?
    photosphere = False             #Use the photosphere as the origin of the photons? If False, the rim of the cocoon will be used

    numericalOpticalDepth = True   #Numerical or analytical optical depth? Analytical integration will assume rho is proportional to R**-s. True: numerical
    analyticRhoProfile = True      #Analytical or manual density profile? True: analytical
    useEATS = True                 #Equal arrival time surface integrator on or off? True: on
#surfaceRings = 15              #If not using EATS, how many rings should the shock front surface be integrated over?
    plotOutput = False              #If running the fitter, and you want to plot the output. Plot: True
    analyseChains = False           #Analyse the output Chains from MulitNest?
    gridStart,gridEnd,gridStep = 12,24,500          #Start of grid (log10), end of grid (log10), and number of grid points
    printProcess = True           #Print the values of the variable values as the fitter is running. Prints when fitter hits a new lowest chi^2
    allowPrint = True             #Allow the programme to return prints? For using in cluster, False is recommended to avoid obscene logs. Crucial messages will be printed if printCrucial == True, no matter this value
    printCrucial = True            #Print crucial messages

#Plot options
    daysOrSec = 'd'                #Plot output in days or seconds? 'd': days,   'h': hours,   'm': minutes   's': seconds.
    fluxAxis = 'Jy'               #Choose flux units for the flux out-put. 'mJy' = milli-Janskys   ;   'mJy' = Janskys

#Dynamics options
    reverseShock = True           #Include a reverse shock component?
    opticalDepth = False           #Take synchrotron self-absorption optical depth into account? True: yes


#Sanity check options
    sanityCheckEATS = True         #Do sanity check of the EATS integrator?
    printStats = False               #Print the characteristics of the jet evolution

#Fit options
    runOption = 'fit'               #Run option. Fitting routine: 'fit'; Plot 3D likelihood surface of two variables: 'surf'; Print lightcurves of values from constants.txt: 'LC'
    chi2_type = 'log'              #Evaluate chi2 logarithmically ('log') or lineraly ('lin')?
    livePoints = 1000               #Number of live points

#Mock observations options


    mockDistribution = 'log'           #Distribution of the mock data points in the temporal grid. For manual time points, enter a float numpy array, e.g. np.array([4.3,7.2e3,4.1e5]). Options: 'log': logarithmic; 'lin': linear
    numberOfMock = 20                 #Number of mock data points. This option is only activated if mockDistribution is set to an automatic option
    mockIntervalTime = [3600,20*86400]##       #Interval of the temporal grid in mock distribution ([lowerValue,upperValue])
    mockIntervalFreq = [1.e10,1.e20]  #See above, for frequency (used when createMock == True and mockDim == 'E' in the runOption == 'LC' mode)
    createMock = False                 #Creates mock observations.
    gaussianError = 0.1               #Gaussian error (1-sigma) in mock distribution
    gaussianOffset = True             #Activates a gaussian offset. This is the only option (for now)
    offsetType = 'log'                #Offset type of mock data. lin - linear offset, log - logarithmic offset
    useData = True                    #True: Uses observations for LC production. False: Uses input data
    plotComponents = True              #Plot thermal, reverse shock and forward shock components? Yes, plot them: True; No, only plot total lightcurve: False
    plotMock = False                   #Plots mock observation. If creating mock observation, this is automatically over-run to TRUE. Recommendation is to keep this FALSE at all times.
    mockDim = 'T'                      #The dimensions to create mock observations in. 'T': time; 'E': energy
    #freqGrid = [5.e14]
#    freqGrid = [5.e14, 2.e18]
    #freqGrid = [4.56e14, 2.4e17]       #R-band and 1 keV
    #freqGrid = [2.4e18]
    #freqGrid = [1.4e9, 4.8e9, 8.4e9, 22.5e9, 4.3e10, 1.e11, 3.e11] #Radio
    #freqGrid = [2.4e18]  #x-rays
    #freqGrid = [ 6.7e14, 8.2e14, 5.4e14, 4.6e14, 3.7e14, 2.5e14, 1.4e14, 4.8e14] #Optical
    #freqGrid = [6.7e14, 8.2e14, 5.4e14, 4.6e14, 3.7e14, 2.5e14, 1.4e14, 4.8e14,1.e16,2.4e18] #xrays-UV-optical
    #freqGrid = [1.e16, 2.4e18] #UV xrays
    #freqGrid = np.concatenate([np.linspace(8.2e14,8.7e14,10) , np.linspace(2.4e17,2.4e18,10)])
    freqGrid = [1.4e9, 4.8e9, 8.4e9, 22.5e9, 4.3e10, 1.e11, 3.e11, 6.7e14, 8.2e14, 5.4e14, 4.6e14, 3.7e14, 2.5e14, 1.4e14, 4.8e14, 2.4e18]   # Will be over-run if (runOption == 'LC') and (createMock == True)
#freqGrid = [1.4e9, 2.4e18]   # Will be over-run if (runOption == 'LC') and (createMock == True)
    tgrid = [1.3e1,1.e8]        #Only used to create mock observations resolved in energies
    
#Parameter constraint options
    fixedRSFSratio = False          #Fixed ratio between the microphysical coefficients of the reverse shock and the forward shock?
    
    
#Output data
    figTypes = 'eps'                                              #File extentions for saving figures
    

#Constants

############################################
#                Intervals                 #
############################################

#outputCounter = 0
    cb = 1.60217657e-12 
    paramLimits = np.array([
            [0.,1.]              #epsilon
            ,[-8.,0.]              #epsilon e (logarithmic)
            ,[-8.,0.]              #epsilon p (logarithmic)
            ,[-8.,0.]              #epsilon e RS (logarithmic)
            ,[-8.,0.]              #epsilon p RS (logarithmic)
            ,[48.,56.]            #log10( E0 ) (isotropic)
            ,[-5.,1.]             #n (logarithmic)
            ,[0.,2.1]            #s
            ,[1.,4.]         #Gamma0 (logarithmic)
            ,[-8.,0.]           #epsilon B (FS, logarithmic)
            ,[-8.,0.]           #epsilon B (RS, logarithmic) 
            ,[1.1,2.8]            #Electron energy distribution slope p for the forward shock
            ,[-5.,5.]             #t0 (logarithmic)
            ,[.0001,90*np.pi/180.]  #theta0 (opening angle)
            ,[0.,20*np.pi/180]            #alpha (off-axis angle)
            ,[0.,3.]              #tN (thermal)
            ,[10.,15.]            #log(r0) (photosphere)
            ,[cb*1.e5,cb*1.e6]    #Tprim0 [eV] (thermal)
            ,[5.,8.]              #log(N0) number of emitted photons
            ,[-1,4]            #Duration of outflow after trigger time (progenitor frame). Only affects the reverse shock scenario. Often assumed to be same as the prompt T_90 time. 10-logarithmic
            ,[1.,3.]        #Electron energy distribution slope p for reverse shock
            ,[8.,14.]           #Initial radius of cocoon
            ,[30.,80.]        #Total number of particles in cocoon (logarithmic)
            ,[1.,5.]            #Duration of outflow producing thermal component (logarithmic)
            ,[0.,80.*np.pi/180]  #Initial opening angle of cocoon
            ,[0,3.]              #Initial Lorentz factor of cocoon
            ,[1e-2,10]     #Redshift
            ,[0.,1.]        #Wv
            ])
##########################################
#        Variables for fitting           #
##########################################

    #Order:       epsil,eps_E,eps_P, eE_RS,eP_RS, E0    n     s   Gamma;epsB; eB_RS   p    t0   theta0 alpha  tN  r0    TMax   N0  prompt, pRS  ,Rcoc N_coc T_coc,thet_C,Gam0_C z    Wv
    parametrar = [False,True ,False,False,False,True ,True ,False,True ,True ,False,True ,False,True ,True ,False,False,False,False,False,False,False,False,False,False,False,False,False] #False = Constant. 


    return inputData,fileExtention,breakName,inputTimeFormat,inputFreqUnit,plotInput,GRBlabel,GRBinput,allowDataRemoval,gridReducer,surfacePlot,thermalComp,photosphere,numericalOpticalDepth,analyticRhoProfile,useEATS,plotOutput,gridStart,gridEnd,gridStep,printProcess,daysOrSec,fluxAxis,reverseShock,opticalDepth,sanityCheckEATS,printStats,runOption,mockDistribution,numberOfMock,mockIntervalTime,mockIntervalFreq,createMock,gaussianError,gaussianOffset,offsetType,useData,plotComponents,plotMock,mockDim,fixedRSFSratio,figTypes,freqGrid,tgrid,cb,paramLimits,parametrar,allowPrint,printCrucial,analyseChains,livePoints,chi2_type
