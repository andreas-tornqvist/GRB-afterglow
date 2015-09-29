def modelFunc(R,constants,fixedRSFSratio,reverseShock,runOption,useEATS,thermalComp,photosphere,printStats,allowPrint,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,createMock,opticalDepth,plotComponents,daysOrSec,fluxAxis,numberOfPoints,printProcess,surfaceRingsIn,ndims,FS_only,RS_only,thermal_only,chi2_type):
    import numpy as np
    from dynamics import dyn
#    from synchEATS import specGenEATS
    from matplotlib import pyplot as plt
    from cosmocalc import cosmocalc
    if runOption=='LC': import time
    from EATS_func import eats_function
    import warnings


    #Constants
    c  ,  pi = 2.9979e10,  np.pi
    mp = 1.6726e-24
    me = 9.1094e-28
    hcgs = 6.6260755e-27   #Planck constant in cgs
    hsi = 6.6260755e-34    #Planck constant in si
    kB = 1.380658e-16
    Bconst = 16.*pi*mp
    mockDim = 'T'
    sigmaT = 6.6524e-25
    qe = 4.803204e-10
    sigma_B = 5.6704e-5  ### Stephan Boltzann constant in cgs
    rad_const = 4 * sigma_B / c   #### Radiation constant
    #kappa13 = -kappa1/3
    #kappa12 = kappa1/2
    #kappa11 = -1/kappa1
    #kappa2p = kappa2*(p-1)/2
    #kappa12inv = -1/kappa2
    #kappa33 = -kappa3/3
    #kappa3p = kappa3*(p-1)/2
    #kappa13inv = -1/kappa3
    #kappa42 = kappa4/2
    #kappa14 = -1/kappa4                                                                                

    epsilon,epsiloneLog,epsilonpLog,epsilone3Log,epsilonp3Log,E0log,nLog,s,Gamma0log,eBlog,eB3log,p,logt0,theta0,alpha,tN,logrphoto,Tprim0,N0,tprompt_log,pRS,R0_Clog,N_Clog,t_outflow_log,theta0_C,GammaC0log,z,WV = constants

    if (theta0 + alpha) < 0.053: surfaceRings = surfaceRingsIn   #If alpha + theta0 < 3 degrees, the number of rings in EATS integrator is constant
    else: surfaceRings = int(surfaceRingsIn * (theta0+alpha) / 0.053) #If alpha + theta0 > 3 degrees, the number of rings in EATS integrator increases linearly



    D = cosmocalc(z,WV=0.0)['DL_cm']    #Open universe: WV = 0, flat universe: WV = 0.73
    #print tprompt

    epsilons_unity = True
    eB , eB3 = 10**eBlog  ,  10**eB3log 
    epsilone = 10**epsiloneLog
    epsilone3 , epsilonp3 = 10**epsilone3Log , 10**epsilonp3Log
    if epsilons_unity: # Energy partitions add up to 1
        epsilonp = 1 - epsilone - eB
    else:
        epsilonp = 10**epsilonpLog
    n = 10 ** nLog
    t0 = -10**logt0
    rphoto = 10**logrphoto
    tprompt = 10**tprompt_log
    Gamma0 = 10**Gamma0log
    
    firstTimePoint = 100 * (1+z) * Gamma0   #The comoving time at 100s observer time
    #tprompt = 20.
    if fixedRSFSratio:
        eB3 = RSFSratio * eB
        epsilone3 = RSFSratio * epsilone
        epsilonp3 = RSFSratio * epsilonp

    
    phipF = 1.89 - 0.935*p + 0.17*p**2
    phipS = 0.54 + 0.08*p
    XpF = 0.455 + 0.08*p
    XpS = 0.06 + 0.28*p
    kappa1 = 2.37 - 0.3*p
    kappa2 = 14.7 - 8.68*p + 1.4*p**2
    kappa3 = 6.94 - 3.844*p + 0.62*p**2
    kappa4 = 3.5 - 0.2*p
   
    #kappas = [-kappa1/3,kappa1/2,-1/kappa1,kappa2*(p-1)/2,-1/kappa2,-kappa3/3,kappa3*(p-1)/2,-1/kappa3,kappa4/2,-1/kappa4,phipF,phipS]
    kappa13 = -kappa1/3
    kappa12 = kappa1/2
    kappa11 = -1/kappa1
    kappa2p = kappa2*(p-1)/2
    kappa12inv = -1/kappa2
    kappa33 = -kappa3/3
    kappa3p = kappa3*(p-1)/2
    kappa13inv = -1/kappa3
    kappa42 = kappa4/2
    kappa14 = -1/kappa4

    #Constants to optical depth calculation

    if opticalDepth:
        alpha0Ffactor = 11.7 * phipF * XpF**(-3) * qe / mp
        alpha0Sfactor = 7.8 * phipS * XpS**(-(4+p)/2.) * (p+2)*(p-1)* qe / mp / (p+2/3.)

   

    if reverseShock:
        phipFRS = 1.89 - 0.935*pRS + 0.17*pRS**2
        phipSRS = 0.54 + 0.08*pRS
        XpFRS = 0.455 + 0.08*pRS
        XpSRS = 0.06 + 0.28*pRS
        kappa1RS = 2.37 - 0.3*pRS
        kappa2RS = 14.7 - 8.68*pRS + 1.4*pRS**2
        kappa3RS = 6.94 - 3.844*pRS + 0.62*pRS**2
        kappa4RS = 3.5 - 0.2*pRS
        
        kappa13RS = -kappa1RS/3
        kappa12RS = kappa1RS/2
        kappa11RS = -1/kappa1RS
        kappa2pRS = kappa2RS*(pRS-1)/2
        kappa12invRS = -1/kappa2RS
        kappa33RS = -kappa3RS/3
        kappa3pRS = kappa3RS*(pRS-1)/2
        kappa13invRS = -1/kappa3RS
        kappa42RS = kappa4RS/2
        kappa14RS = -1/kappa4RS

        if opticalDepth:
            alpha0FRSfactor = 11.7 * phipFRS * XpFRS**(-3) * qe / mp
            alpha0SRSfactor = 7.8 * phipSRS * XpSRS**(-(4+pRS)/2.) * (pRS+2)*(pRS-1)* qe / mp / (pRS+2/3.)

    profile_cutoff = 18

    E0 = 10**E0log * (1 - np.cos(theta0)) / 2 #E0log is the isotropic energy, here it is corrected for the opening angle theta0
    M0 = E0*c**-2 / Gamma0
    M03 = 0.#2*pi*R[0]**3*(1-np.cos(theta0)) * n * mp /3              #M0/1000
    twoPi = 2*pi
    A0 = n * (mp + me) * 10**(profile_cutoff*s)
    Rdec = ((3-s)*M0/(4*pi*A0*Gamma0)) ** (1/(3.-s))
    mmean = (me+mp)/2.   #Constant
    

    
    #Get lightcurve from model. Exctract points corresponding to data. Compare (get chi**2) and return the log-likelihood
    
    #Find the last time that we want to calculate model for. 
    tobsEnd = np.max(tdata + t0)
    if runOption=='fit': tobsRedUpper = tobsEnd
    else:
        tobsRedLower = 1.
        tobsRedUpper = 2.2e10 #250 days
        tobsEnd = np.copy(tobsRedUpper)

    if (runOption=='LC'):
        dynStartTime = time.time()

    if reverseShock: Gamma,tburst,tobs,tcomoving,theta,M,m,M3,rho4,RScutOff,Rcut,gamma43,gammae,gammaeRS = dyn(R,Gamma0,epsilone,epsilonp,epsilone3,epsilonp3,eB,eB3,E0,M0,M03,s,A0,t0,p,pRS,theta0,z,tprompt,reverseShock,tobsRedUpper*2,profile_cutoff)
    else: Gamma,tburst,tobs,tcomoving,theta,M,m,Rcut,gammae = dyn(R,Gamma0,epsilone,epsilonp,epsilone3,epsilonp3,eB,eB3,E0,M0,M03,s,A0,t0,p,pRS,theta0,z,tprompt,reverseShock,tobsRedUpper*2,profile_cutoff)
    
    
    if False:  #Plotting the lorentz factor
        gammaWrite = open('gamma.txt','w')
        gammaWrite.write(' '.join(map(str,Gamma)))
        gammaWrite.close()
        tobsWrite = open('tobs.txt','w')
        tobsWrite.write(' '.join(map(str,tobs)))
        tobsWrite.close()
        


    if thermalComp:
        cocSurfaceRings = 20
        tau_photo = 1.
        r0_coc = 1e7   ### Initial guess. This is the radius of the jet origin (basically the Schwarzschild radius of the star) in cm
        kappa_0 =  5e24
        ### Allocating space
        R_photo_int = np.zeros(cocSurfaceRings)
        T_int = np.zeros(cocSurfaceRings)

        precalc_length = 100   ### Length of the grid used for precalculation of cocoon outflow properties
        ### See note from 20/6 2015
        
        t_outflow = 10**t_outflow_log  ### The duration of the outflow
        t_start = np.linspace(0,t_outflow,precalc_length)  ### The grid of the starting time (prog. time) of all outflow segments
        R0_coc = 10**R0_Clog ### Starting radius of cocoon
        MC = 10**N_Clog * mp
#        RcocIn = np.logspace(R0_Clog,R0_Clog+6,1000)
#        R0_C = 10**R0_Clog
        GammaC = 10**GammaC0log#np.sqrt(3/2.)
        beta_coc = np.sqrt(1-1/GammaC**2)


        ### Constructing temporal grids
        t0_coc = R0_coc / beta_coc / c
        last_t_co_coc = tobsEnd * (1+beta_coc) * GammaC / (1+z) + t0   ### The comoving time when the last outflow segment passes tobsEnd
        t_co_rel = np.logspace(np.log10(t0_coc),np.log10(last_t_co_coc-t0),precalc_length)  ###The relative elapsed comoving time since outflow passed radius R0_coc
    
        #tobs_coc_rel = (1+z) / (1+beta_coc)/GammaC * t_co_rel  ### The relative observer time elapsed since outflow passed radius R0_coc

        ### Radial grid
        R_coc = np.logspace(np.log10(R0_coc),np.log10(R0_coc + last_t_co_coc * beta_coc * c),precalc_length) ### Constructed from burster frame

        ### Estimation of opening angle
        gammabar_coc = (4*GammaC+1) / (3*GammaC)
        vperp_coc = ( gammabar_coc*(gammabar_coc-1)*(GammaC-1) / (1 + gammabar_coc*(GammaC-1)) )**(1/2.) * c
        theta_C = np.zeros(precalc_length)
        theta_C[0] = theta0_C
        for i_theta_coc in range(1,precalc_length):
            theta_C[i_theta_coc] = theta_C[i_theta_coc-1] + vperp_coc / R_coc[i_theta_coc] * (t_co_rel[i_theta_coc] - t_co_rel[i_theta_coc-1])
            if theta_C[i_theta_coc] >= np.pi/2:
                theta_C[i_theta_coc:] = np.pi/2
                break

        fovAngleCoc = 1 / GammaC
        fovInRimCoc = fovAngleCoc<(theta_C+alpha)
        fovInRimCoc_lower = fovAngleCoc<(theta_C - alpha)
        lowerLimitCoc = (theta_C - alpha) * (fovInRimCoc_lower==False) + fovAngleCoc * fovInRimCoc
        upperLimitCoc = (alpha+theta_C) * (fovInRimCoc==False) + fovAngleCoc * (fovInRimCoc)   ### The upper angle that is visible to the observer
        #tobsRimC = (1+z) * (tburstC - R_coc * np.cos(fovAngleCoc * (fovInRimCoc)  +  (theta_C+alpha) * (fovInRimCoc==False)) / c)     # Observing time at the rim for each radial point. Will use this in setting EATS grid. Note from 9/12 -13


        T = (E0 / 4 / np.pi / r0_coc**2 / c / rad_const)**(1/4.) / GammaC * (R_coc/r0_coc)**(-2/3)   ### Precalculation. After Pe'er 2007.
        rho_coc = MC / precalc_length / (2*np.pi * (1-np.cos(theta_C)) * R_coc**2)   ### Precalculation of density
        tau_spec = kappa_0 * rho_coc**2 * T**(-7/2.)   ### Precalculation of the specific optical depth (tau/delta R)
        tau_segment = tau_spec[1:] * (R_coc[1:] - R_coc[:-1])  ### The optical depth the shortest way in each radial segment. Used when finding photosphere in EATS integrator
        ### Creating a grid of location of the photosphere on the viewing angle axis. Moving away from the viewing angle, the photosphere will move to a larger radius R_coc
        R_photo_pre = np.zeros([precalc_length,cocSurfaceRings])
        R_photo_pre[0] = np.copy(R_coc[0])
        for i_photo_coc in range(1,precalc_length):
            ring_segments_coc = np.linspace(0,upperLimitCoc[i_photo_coc],cocSurfaceRings)
            for i_phococ_ringseg in range(cocSurfaceRings):
                tau_pre = 0.
                for i_ph_int in np.linspace(i_photo_coc,1,i_photo_coc).astype(int):
                    try:tau_pre += tau_segment[i_ph_int-1] / np.cos(ring_segments_coc[i_phococ_ringseg])   ### Integrating optical depth
                    except: 
                        print tau_pre
                        print i_ph_int
                    if tau_pre >= tau_photo:   ### Has the optical depth reached the value where gas is no longer transparent?
                        R_photo_pre[i_photo_coc,i_phococ_ringseg] = np.copy(R_coc[i_ph_int-1])
                        break
            if tau_pre < tau_photo: R_photo_pre[i_photo_coc] = np.copy(R_coc[0])
        tobs_conv_fac = (1-beta_coc) / beta_coc / c
        tobs_coc_photo = (t_start + (R_photo_pre[:,0]-R0_coc)*tobs_conv_fac) * (1+z)   ### The observer time to each location of the photosphere along the line of sight

#        plt.figure(1)
#        plt.plot(R_coc,R_photo_pre)
#        plt.loglog()
#        plt.figure(2)
#        plt.plot(R_coc,rho_coc)
#        plt.loglog()
#        plt.figure(3)
#        plt.plot(R_coc,tau_spec)
#        plt.loglog()
#        plt.figure(4)
#        plt.plot(R_coc,theta_C)
#        plt.xscale('log')
#        plt.figure(5)
#        plt.plot(R_coc,tobs_coc_photo)
#        plt.loglog()
#        plt.figure(6)
#        plt.plot(R_coc,t_start)
#        plt.loglog()
#        plt.show()


#        T0 = None
#        GammaC,tburstC,tobsC,tcomovingC,theta_C,_,_,_,_,Tcoc,Rcoc = dyn(RcocIn,GammaC0,T0,epsilonp,epsilone3,epsilonp3,eB,eB3,E0,MC,M03,s,0.,t0,theta0_C,z,tprompt,False,tobsEnd,profile_cutoff)
#        if True: ### Setting the initial temperature through the relation of Pe'er et al. (2007)
#            T = (E0 / 4 / np.pi / r0**2 / c / rad_const) ** 0.25 / GammaC * (Rcoc / r0) ** (-2/3.)
#        else:
#            T0 = 10**Tlog
#        print Tcoc
#        plt.plot(Rcoc,Tcoc)
#        plt.loglog()
#        plt.show()
#        cocSurfaceRings = 500



#        try:
#            if not np.isnan(GammaC)[-1]:
#                cocCrash = False
#                beta_coc = np.sqrt(1-1/GammaC**2)
#                angCocInd = np.zeros(cocSurfaceRings,dtype=int)
#            else: cocCrash = True
#        except: cocCrash = True
#        while cocCrash == True:
#            print 'Cocoon dynamics module crashed. Tightening grid.'
#            try: tighterCoc *= 2
#            except: 
#                tighterCoc = 1000
#                cocIterator = 0
#            GammaC,tburstC,tobsC,tcomovingC,theta_C,_,_,_,_,Tcoc,Rcoc = dyn(RcocIn,GammaC0,T0,epsilonp,epsilone3,epsilonp3,eB,eB3,E0,MC,M03,s,0.,t0,theta0_C,z,tprompt,False,tobsEnd,profile_cutoff)
#            try: 
#                if not np.isnan(GammaC)[-1]:   #If cocoon dynacism calculation worked without crashing
#                    cocCrash = False
#                    beta_coc = np.sqrt(1-1/GammaC**2)
#                    angCocInd = np.zeros(cocSurfaceRings,dtype=int)
#            except:
#                RcocIn = np.logspace(R0_Clog,R0_Clog+6,1000*tighterCoc)
#                cocIterator += 1
#                if cocIterator >= 1: 
#                    thermalComp = False
                    
#                    break
#        if cocCrash == True: print 'Cocoon crashed. Running this point without thermal'
#        else:
#            if len(Rcoc) != len(GammaC):
#                Rcoc = Rcoc[:len(GammaC)]
#            where_kill_cocoon = np.where(Tcoc<1e4)#np.where(Dtau[1:]/(Rcoc[1:]-Rcoc[0]) >= 1)
#            try:
#                kill_cocoon = np.min(where_kill_cocoon)+1
#            except:
#                kill_cocoon = len(tobsC)-1
#            if runOption=='LC':print 'Cocoon becomming transparent at %d s (observer time)'%tobsC[kill_cocoon]
#            if photosphere:
#                Rcoc_dist = np.zeros(len(Rcoc))
#                Rcoc_dist[0] = (Rcoc[1] - Rcoc[0]) / 2
#                Rcoc_dist[-1] = (Rcoc[-1] - Rcoc[-2]) / 2
#                Rcoc_dist[1:-1] = (Rcoc[2:] - Rcoc[:-2]) / 2
#                rho_coc = np.zeros(len(Rcoc))
#                for i_coc in range(1,len(rho_coc)):   ### Calculating the density of the cocoon for each radial step. Only used when photosphere==True
#                    rho_coc[i_coc] = MC / (2 * np.pi * np.sum((Rcoc[1:i_coc+1]-Rcoc[0])**2*Rcoc_dist[1:i_coc+1]*(1-np.cos(theta_C[1:i_coc+1]))))  
                #print 2 * np.pi * np.sum((Rcoc[:i_coc+1]-Rcoc[0])**2*Rcoc_dist[:i_coc+1]*(1-np.cos(theta_C[:i_coc+1])))
                #print 
                #raw_input(rho_coc[i_coc])
#                rho_coc[0] = np.copy(rho_coc[1])
                
                #plt.show()
#                tau_coc = 1  ### At which optical depth the photosphere exists
                #Setting these two first lines to see evaluation
#                kappa_0 =  5e24
#                kappa_R = kappa_0 * rho_coc * Tcoc**(-7/2.)
                
                #plt.show()
#                Dtau = tau_coc / rho_coc / kappa_R   ### The distance with optical depth tau=1
                #plt.plot(tobsC,Dtau)
                #plt.loglog()
                #plt.show()
                #plt.plot(tburstC,Tcoc)
                #plt.loglog()
                #plt.show()
                ### Find where the cocoon becomes transparent
                
            

    try: 
        while np.isnan(Gamma): #Routine to tighten the grid if the dynamics module crashed
            print len(R)
            try: crashedOnce += 1
            except: crashedOnce = 1
            R = np.logspace(np.log10(R[0]),np.log10(R[-1]),len(R)*2)
            if reverseShock: Gamma,tburst,tobs,tcomoving,theta,M,m,M3,rho4,RScutOff,Rcut,gamma43,gammae,gammaeRS = dyn(R,Gamma0,epsilone,epsilonp,epsilone3,epsilonp3,eB,eB3,E0,M0,M03,s,A0,t0,p,pRS,theta0,z,tprompt,reverseShock,tobsRedUpper*2,profile_cutoff)
            else: Gamma,tburst,tobs,tcomoving,theta,M,m,Rcut,gammae = dyn(R,Gamma0,epsilone,epsilonp,epsilone3,epsilonp3,eB,eB3,E0,M0,M03,s,A0,t0,p,pRS,theta0,z,tprompt,reverseShock,tobsRedUpper*2,profile_cutoff)

            if crashedOnce >= 6: break
    except: 
        if (runOption=='LC'): print "Dynamics module time use: %f"%(time.time()-dynStartTime)
    try: len(Gamma)
    except:
#        for print_crash in range(n_param):
        print 'density = %s'%n
        print 'Energy = %s'%E0
        print 'Gamma0 = %s'%Gamma0
            
        print "\nDynamics module crashed once!\n"
        return float('nan'),[],[]
    #Rcut = R[:len(Gamma)]        #Cleaning up, since dynamics module stop when tobs > tobsEnd
    gammaAd = (4 + 1/Gamma) / 3
    GammaEff = (gammaAd * Gamma**2 - gammaAd + 1) / Gamma
    
    
    beta = np.sqrt(1-1/Gamma**2)        
    rho = A0*Rcut**-s                                                                                                                               
    nprim = 4*rho*Gamma/mp
    
#    plt.plot(Rcut,Gamma*beta)
#    plt.loglog()
#    plt.xlabel('tobs')
#    plt.ylabel(r'$\Gamma\beta$')
#    plt.show()

    if False: #Trial of a grid reducer
        reduceSize = 50
        Rcut,Gamma,tburst,tobs,tcomoving,theta,M,M3,rho4 =         R[map(int,np.linspace(1,len(Gamma)-1,reduceSize))],Gamma[map(int,np.linspace(1,len(Gamma)-1,reduceSize))],tburst[map(int,np.linspace(1,len(Gamma)-1,reduceSize))],tobs[map(int,np.linspace(1,len(Gamma)-1,reduceSize))],tcomoving[map(int,np.linspace(1,len(Gamma)-1,reduceSize))],theta[map(int,np.linspace(1,len(Gamma)-1,reduceSize))],M[map(int,np.linspace(1,len(Gamma)-1,reduceSize))],M3[map(int,np.linspace(1,len(Gamma)-1,reduceSize))],rho4[map(int,np.linspace(1,len(Gamma)-1,reduceSize))]


    B = (Bconst*eB*nprim*Gamma*(gammae-1))**(1/2.)*c
    PmaxF = phipF*2.234*qe**3*nprim*B/me/c**2
    PmaxS = phipS*11.17*(p-1)*qe**3*nprim*B/(3*p-1)/me/(c**2)



    

    
    cosTheta = np.cos(theta)

    #Injected material density in progenitor frame
   # tpromptProg = (tprompt-t0)*(1+z)
   # delta0prog = c*tpromptProg
   # rho4prog = M0 / (4*pi*R**2*delta0prog)

    Nterm = (twoPi*Rcut**2*(1-cosTheta))
    N = M / Nterm  / mmean   #Column number density
    fovAngle = 1 / Gamma    #Field of view angle. Define better!
    fovInRim = fovAngle<(theta+alpha)
    chi2 = 0.
    timeGrid = 100   #How tight the lightcurve grid should be when runOption=='LC'
    timeGridSigma = 50



    if (runOption == 'LC'):
        if createMock:
            lightcurve = np.zeros(np.shape(tdata))
            tobsGrid = []
        else:
            lightcurve = np.zeros([iterationLength,timeGrid])
            tobsGrid = np.zeros([iterationLength,timeGrid])

    if runOption == 'one-sigma': #Evaluating lightcurves to plot one-sigma
        lightcurve = np.zeros([iterationLength,timeGridSigma])
        tobsGrid = np.zeros([iterationLength,timeGridSigma])


    if (runOption=='LC'): dynTimeStart = time.time()

#### Spectrum generation ####

    
    ## Setting spectrum generation constants and vectors
    
    gammamFac1 = (p-2)/(p-1)
    gammamFac2 = epsilone*mp/me
    gammam = gammamFac1 * (gammamFac2 * (gammae -1) + 1)
    
    pimec = 2 * pi * me * c
    gammac1t = 3*pimec/sigmaT/tobs
    gammac = gammac1t/Gamma/B**2
    gamToNuFactor = qe * B / pimec
    num = gammam**2 * gamToNuFactor
    nuc = gammac**2 * gamToNuFactor


    twoPiArray = np.array([twoPi]*surfaceRings)
    PprimTemp = np.zeros(surfaceRings)
    phiInter = np.zeros(surfaceRings)
    angleInd = np.zeros(surfaceRings,dtype=int)
    
    if reverseShock: 
        RS_elements = np.count_nonzero(rho4)
        gamma43Ad = (4 + 1/gamma43[:RS_elements])
        gamma43Eff = (gamma43Ad * gamma43[:RS_elements]**2 - gamma43Ad + 1) / gamma43[:RS_elements]
        rho3prim = 4*rho4[:RS_elements]*gamma43[:RS_elements]
        n3prim = rho3prim / mmean
        BRS = (Bconst*eB3*n3prim*Gamma[:RS_elements]*(gammae[:RS_elements]-1))**(1/2.)*c
        PmaxF_RS = phipF*2.234*qe**3*n3prim*BRS/me/c**2
        PmaxS_RS = phipS*11.17*(p-1)*qe**3*n3prim*BRS/(3*pRS-1)/me/(c**2)

        gammamFac1RS = (pRS-2)/(pRS-1)
        gammamFac2RS = epsilone3*mp/me
        gamToNuFactorRS = qe * BRS / pimec
        gammamRS = gammamFac1RS * (gammamFac2RS * (gammaeRS[:RS_elements] -1) + 1)

        gammacRS = gammac1t[:RS_elements]/ gamma43[:RS_elements] / BRS**2
#        numRS = -np.ones(len(Gamma))
#        nucRS = -np.ones(len(Gamma))
        numRS = gammamRS**2 * gamToNuFactorRS
        nucRS = gammacRS**2 * gamToNuFactorRS



#    Pforward , Pbehind = np.zeros(surfaceRings) , np.zeros(surfaceRings)

    tobsRim = (1+z) * (tburst - Rcut * np.cos(fovAngle * (fovInRim)  +  (theta+alpha) * (fovInRim==False)) / c)     # Observing time at the rim for each radial point. Will use this in setting EATS grid. Note from 9/12 -13

    for nuIte in range(iterationLength):     #Loop over all input frequencies or times

        freqArr = np.array([freq[nuIte]])
        tobsRed = tdata[nuIte,:numberOfEmpties[nuIte]] + t0
        noChi2 = ((runOption == 'LC') & (createMock)) | (runOption == 'one-sigma')  #This is true if we want to produce mock observations, and don't want to read in data
        if not noChi2:
            if (runOption=='LC') and (not createMock):  #Creating an equally spaced temporal grid to make smoother plots
                tdataLC = tobsRed
                tobsRed = np.logspace(np.log10(tobsRedLower),np.log10(tobsRedUpper),timeGrid)
                tobsGrid[nuIte] = tobsRed
            Fdata = FdataInput[nuIte,:numberOfEmpties[nuIte]]
            errorbar = errorbarInput[nuIte,:numberOfEmpties[nuIte]]
        elif runOption == 'one-sigma': #If plotting the one-sigma range
            tobsRed = np.logspace(np.log10(tobsRedLower),np.log10(tobsRedUpper),timeGridSigma)
            tobsGrid[nuIte] = tobsRed
        if useEATS:
            
        #Equal Arrival Time Surface (EATS) integrator
            
            

            #Allocating space
            
#            if mockDim == 'T':
            
            EATSsteps = len(tobsRed)
#                freqNow = freqArr
#            elif mockDim == 'E':  #Option 'E' is obsolete
#                EATSsteps = len(freqArr)
#                tobsRimNow = tobsRed
    
            Pprim = np.zeros(EATSsteps)     
            if reverseShock: reverseComponent = np.zeros(EATSsteps)
            if thermalComp: 
                thermal_component = np.zeros(EATSsteps)

            for rimI in range(EATSsteps):  #If we want a frequency resolution, len(tobsRed) = 1 . Note that len(tobsRed) must be a vector still, but with one element. Test: Run a file with only one point! 
                #Finding index for the nearest point behind seaked rim radius.

                tobsRimNow = tobsRed[rimI]     #Tells the program what observer's time we want to find emission for
                indRim = np.argmin(np.abs(tobsRimNow-tobsRim))    #Find what index the point at the rim has with the observer's time we are looking for.
                indRim -= (tobsRim[indRim] > tobsRimNow)          #Making sure found index is behind seaked point
                
                
                if tobsRimNow==0:raw_input(tobsRimNow)
                if tobs[indRim] == 0: raw_input(tobs[indRim])
                if tobs[indRim+1] == 0: raw_input(tobs[indRim+1])
                    #Weights for the rim interpolation
                weightRim,weight1Rim,weight2Rim = np.log10(tobs[indRim+1] / tobs[indRim])   ,   np.log10(tobs[indRim+1] / tobsRimNow)   ,   np.log10(tobsRimNow / tobs[indRim])
                if np.isnan(weightRim):
                    print tobs
                    print Gamma[indRim]
                    print Gamma[indRim+1]
                    print tobs[indRim+1]
                    print tobs[indRim]
                    print indRim
                    print len(tobs)
                    raw_input('1')
                if np.isnan(weight1Rim):raw_input('2')
                if np.isnan(weight2Rim):raw_input('3')
                    #Interpolating all surface rings
                if fovInRim[indRim]:  ### Field of view is inside the rim of the jet
                    thetaRimPre = np.copy(fovAngle[indRim])
                else:
                    thetaRimPre = np.copy(theta[indRim])
                
                if fovInRim[indRim+1]:  ### Field of view is inside the rim of the jet
                    thetaRimPre2 = np.copy(fovAngle[indRim+1])
                else:
                    thetaRimPre2 = np.copy(theta[indRim+1])


#                thetaRimPre = theta[indRim] * (fovInRim[indRim]==False) + fovAngle[indRim] * (fovInRim)  #The field of view angle (if the fov covers the jet, fov=theta)
#                thetaRimPre2 = theta[indRim+1] * (fovInRim[indRim+1]==False) + fovAngle[indRim+1] * (fovInRim[indRim+1])
                thetaRim = (thetaRimPre * weight1Rim + thetaRimPre2 * weight2Rim) / weightRim
                lowerLimit = thetaRim - alpha
                upperLimit = thetaRim + alpha

                #Allocating space and setting variable angular step size
                thetaRimDeg = thetaRim * 180 / pi
		#Approximating number of rings in the EATS integration

                



                #GammaInter,mInter,BInter,RInter,tcoInter,thetaInter,phiInter,rho4progInter,tobsInter = np.zeros(surfaceRings),np.zeros(surfaceRings),np.zeros(surfaceRings),np.zeros(surfaceRings),np.zeros(surfaceRings),np.zeros(surfaceRings),np.zeros(surfaceRings),np.zeros(surfaceRings),np.zeros(surfaceRings)
                #if reverseShock: M3Inter,BRSInter,gamma43Inter = np.zeros(surfaceRings),np.zeros(surfaceRings),np.zeros(surfaceRings)

                



                    #Angular grid
                xAngInt = np.linspace(0,thetaRim+alpha,surfaceRings+1)    #The shell is divided into surfaceRing number of rings, with surfaceRings+1 number of borders. xAnd has the central angle of each ring segment, while xAngInt has the border angle. 
                xAng = (xAngInt[:-1] + xAngInt[1:]) / 2
                

                ### Use interpolating when calculating forward shock?
                F_interpolation = True

                
                cosAng = np.cos(xAng)
                 ### Looping along the Equal Arrival Time Surface

                try:angleInd[0] = np.argmin(np.abs(tobs-tobsRed[rimI]))  #Index of the time at the surface point in the observer's LoS
                except:
                    print surfaceRings
                    print 'Oooops angleInd[0] failed!!!'
                    print surfaceRingsIn
                    print theta0
                    print alpha
                    print len(tobsRed)
                    print len(tobs)
                    print len(angleInd)



                for indAng in range(1,surfaceRings):
                    angleInd[indAng] = np.copy(angleInd[indAng-1])
                    while tobsRed[rimI] < ((1+z) * (tburst[angleInd[indAng]] - Rcut[angleInd[indAng]] * cosAng[indAng] / c)):
                        if angleInd[indAng] <= 0:  ### If a too early point is seeked
                            break
                        angleInd[indAng] -= 1  #Finding the index corresponding to the shell behind the currently evaluated point on the EATS


                if reverseShock:
                    where_angleInd_RS = np.where(angleInd < (RS_elements-1))
                    angleInd_RS = np.copy(angleInd[where_angleInd_RS])   ### These elements should be used when calculating radiation from RS


                tobsBehind = (1+z) * (tburst[angleInd] - Rcut[angleInd] * cosAng / c)
                phiInter = (xAng < lowerLimit) * twoPiArray

                if alpha != 0: #Off-axis
                    partialRingsInd = np.where((xAng >= lowerLimit) & (xAng < upperLimit))   #Rings crossing the rim. Happens only when alpha != 0
                    offAxisFoV = (theta[angleInd[partialRingsInd]]**2 - alpha**2 - xAng[partialRingsInd]**2) / (2*alpha*xAng[partialRingsInd])
                    offAxisFoV[np.where(offAxisFoV<-1)] = -1.
                    if np.sum(offAxisFoV > 1) != 0: offAxisFoV[np.where(offAxisFoV > 1)] = 1.
                    phiInter[partialRingsInd] = 2*np.pi - 2*np.arccos(offAxisFoV)


                onePzFreq = (1+z) * freqArr
                cosAng2 = np.cos(xAng)**2
                if F_interpolation:
                    nuPrimForward = onePzFreq * Gamma[angleInd+1] * (1-beta[angleInd+1] * np.cos(xAng))
                    nuPrimBehind = onePzFreq * Gamma[angleInd] * (1-beta[angleInd] * np.cos(xAng))

                    rhoPrimBehind = 4 * rho[angleInd] * Gamma[angleInd]
                    rhoPrimForward = 4 * rho[angleInd+1] * Gamma[angleInd+1]

                    RForward2 = Rcut[angleInd+1]**2
                    RBehind2 = Rcut[angleInd]**2
                                                                                                          
                    thickForward = m[angleInd+1] / (twoPi*(1-np.cos(theta[angleInd+1])) * Gamma[angleInd+1] * rhoPrimForward * RForward2)
                    thickBehind = m[angleInd] / (twoPi*(1-np.cos(theta[angleInd])) * Gamma[angleInd] * rhoPrimBehind * RBehind2)
                
                    
                

                #Forward shock

                    slow_cooling_low = np.where(num[angleInd] <= nuc[angleInd])
                    slow_cooling_high = np.where(num[angleInd+1] <= nuc[angleInd+1])

                    fast_cooling_low = np.where(num[angleInd] > nuc[angleInd])
                    fast_cooling_high = np.where(num[angleInd+1] > nuc[angleInd+1])

                    Pbehind = eats_function(p, nuPrimBehind, slow_cooling_low, fast_cooling_low, num[angleInd], nuc[angleInd], PmaxF[angleInd], PmaxS[angleInd])
                    Pforward = eats_function(p, nuPrimForward, slow_cooling_high, fast_cooling_high, num[angleInd+1], nuc[angleInd+1], PmaxF[angleInd+1], PmaxS[angleInd+1])

                    Pforward *= 1e23 * RForward2 * thickForward / (Gamma[angleInd+1] ** 2 * (1-beta[angleInd+1]*cosAng)**2)
                    Pbehind *=  1e23 * RBehind2 * thickBehind / (Gamma[angleInd] ** 2 * (1-beta[angleInd]*cosAng)**2)

                    tobsFront = (1+z) * (tburst[angleInd+1] - Rcut[angleInd+1] * np.cos(xAng) / c)
                        #Weights for the intermediate angle interpolations. Tested, works fine
                    weight,weight1,weight2 = np.log10( tobsFront / tobsBehind )  ,  np.log10( tobsFront / tobsRimNow )  ,  np.log10( tobsRimNow / tobsBehind )

                    PprimTemp = (Pbehind * weight1 + Pforward * weight2)  /  weight  
                    #Reverse shock
#                    if reverseShock: 

#                        RForward2RS = Rcut[angleInd_RS+1]**2
#                        RBehind2RS = Rcut[angleInd_RS]**2

#                        cosAngRS = np.cos(xAng[where_angleInd_RS])
#                    cosAng2 = np.cos(xAng[np.where(rho4[angleInd]!=0)])**2

#                        rho3primBehind = rho3prim[angleInd_RS]
#                        rho3primForward = rho3prim[angleInd_RS+1]
#                    PRSforward = np.zeros(surfaceRings)
#                    PRSbehind = np.zeros(surfaceRings)

#                        slow_cooling_lowRS = np.where(numRS[angleInd_RS] <= nucRS[angleInd_RS])
#                        slow_cooling_highRS = np.where(numRS[angleInd_RS+1] <= nucRS[angleInd_RS+1])

#                        fast_cooling_lowRS = np.where(numRS[angleInd_RS] > nucRS[angleInd_RS])
#                        fast_cooling_highRS = np.where(numRS[angleInd_RS+1] > nucRS[angleInd_RS+1])

#                        print slow_cooling_lowRS
#                        print slow_cooling_highRS
#                        print fast_cooling_lowRS
#                        print fast_cooling_highRS
#                        print numRS[angleInd_RS+1]
#                        raw_input(numRS[angleInd_RS])


#                        PRSforward = eats_function(pRS, nuPrimForward, slow_cooling_highRS, fast_cooling_highRS, numRS[angleInd_RS+1], nucRS[angleInd_RS+1], PmaxF_RS[angleInd_RS+1], PmaxS_RS[angleInd_RS+1])
#                        PRSbehind = eats_function(pRS, nuPrimBehind, slow_cooling_lowRS, fast_cooling_lowRS, numRS[angleInd_RS], nucRS[angleInd_RS], PmaxF_RS[angleInd_RS], PmaxS_RS[angleInd_RS])


            #        PRSforward = (numRS[angleInd+1] > nucRS[angleInd+1]) * PmaxF_RS[angleInd+1] * ((nuPrimForward/nucRS[angleInd+1])**(kappa13RS) + (nuPrimForward/nucRS[angleInd+1])**(kappa12RS)) ** (kappa11RS) * (1+(nuPrimForward/numRS[angleInd+1])**(kappa2pRS))**(kappa12invRS) + (numRS[angleInd+1] <= nucRS[angleInd+1]) * PmaxS_RS[angleInd+1] * ((nuPrimForward/numRS[angleInd+1])**(kappa33RS) + (nuPrimForward/numRS[angleInd+1])**(kappa3pRS))**(kappa13invRS) * (1+(nuPrimForward/nucRS[angleInd+1])**(kappa42RS))**(kappa14RS)
                    
            #        PRSbehind = (numRS[angleInd] > nucRS[angleInd]) * PmaxF_RS[angleInd] * ((nuPrimBehind/nucRS[angleInd])**(kappa13RS) + (nuPrimBehind/nucRS[angleInd])**(kappa12RS)) ** (kappa11RS) * (1+(nuPrimBehind/numRS[angleInd])**(kappa2pRS))**(kappa12invRS) + (numRS[angleInd] <= nucRS[angleInd]) * PmaxS_RS[angleInd] * ((nuPrimBehind/numRS[angleInd])**(kappa33RS) + (nuPrimBehind/numRS[angleInd])**(kappa3pRS))**(kappa13invRS) * (1+(nuPrimBehind/nucRS[angleInd])**(kappa42RS))**(kappa14RS)

                    

#                        thickRSForward = M3[angleInd_RS+1]/ (4*twoPi*(1.-np.cos(theta[angleInd_RS+1])) * gamma43_full[angleInd_RS+1]**2*rho4[angleInd_RS+1]*RForward2RS)
#                        thickRSBehind = M3[angleInd_RS]/ (4*twoPi*(1.-np.cos(theta[angleInd_RS])) * gamma43_full[angleInd_RS]**2*rho4[angleInd_RS]*RForward2RS)

#                        PRSforward *= 1e23 * RForward2RS * thickRSForward / (Gamma[angleInd_RS+1] ** 2 * (1-beta[angleInd_RS+1]*cosAngRS)**2)
#                        PRSbehind *=  1e23 * RBehind2RS * thickRSBehind / (Gamma[angleInd_RS] ** 2 * (1-beta[angleInd_RS]*cosAngRS)**2)
                
#                        tobsBehindRS = (1+z) * (tburst[angleInd_RS] - Rcut[angleInd_RS] * cosAngRS / c)
                    
#                        tobsFrontRS = (1+z) * (tburst[angleInd_RS+1] - Rcut[angleInd_RS+1] * cosAngRS / c)
                        
#                        weightRS,weight1RS,weight2RS = np.log10( tobsFrontRS / tobsBehindRS )  ,  np.log10( tobsFrontRS / tobsRimNow )  ,  np.log10( tobsRimNow / tobsBehindRS )



#                        PRSprimTemp = (PRSbehind * weight1RS + PRSforward * weight2RS)  /  weightRS

                        
                else:
                    
                    nuPrim = onePzFreq * Gamma[angleInd] * (1-beta[angleInd] * np.cos(xAng))
                    rhoPrim = 4 * rho[angleInd] * Gamma[angleInd]
                    R2 = Rcut[angleInd]**2
                    thickness = m[angleInd] / (twoPi*(1-np.cos(theta[angleInd])) * Gamma[angleInd] * rhoPrim * R2)

                #Forward shock

                    slow_cooling = np.where(num[angleInd] <= nuc[angleInd])
                    fast_cooling = np.where(num[angleInd] > nuc[angleInd])
                    P = eats_function(p, nuPrim, slow_cooling, fast_cooling, num[angleInd], nuc[angleInd], PmaxF[angleInd], PmaxS[angleInd])
                    PprimTemp =  P*1e23 * R2 * thickness / (Gamma[angleInd] ** 2 * (1-beta[angleInd]*cosAng)**2)

                    #Reverse shock
                if reverseShock: 
                    #nuPrim = onePzFreq * Gamma[angleInd] * (1-beta[angleInd] * np.cos(xAng))
#                    nuPrimRS = onePzFreq * Gamma[angleInd] * (1-beta[angleInd] * np.cos(xAng))
                    RBehind2_RS = Rcut[angleInd_RS]**2
                    RForward2_RS = Rcut[angleInd_RS+1]**2

#                    cosAngRS = np.cos(xAngRS)
                    slow_coolingRS_behind = np.where(numRS[angleInd_RS] <= nucRS[angleInd_RS])
                    slow_coolingRS_forward = np.where(numRS[angleInd_RS+1] <= nucRS[angleInd_RS+1])


                    fast_coolingRS_behind = np.where(numRS[angleInd_RS] > nucRS[angleInd_RS])
                    fast_coolingRS_forward = np.where(numRS[angleInd_RS+1] > nucRS[angleInd_RS+1])


                    thicknessRS_behind = M3[angleInd_RS]/ (4*twoPi*(1.-np.cos(theta[angleInd_RS])) * gamma43[angleInd_RS]**2*rho4[angleInd_RS]*RBehind2_RS)
                    thicknessRS_forward = M3[angleInd_RS+1]/ (4*twoPi*(1.-np.cos(theta[angleInd_RS+1])) * gamma43[angleInd_RS+1]**2*rho4[angleInd_RS+1]*RForward2_RS)

                    PRS_b_fac = 1e23 * RBehind2_RS * thicknessRS_behind / (Gamma[angleInd_RS] ** 2 * (1-beta[angleInd_RS]*cosAng[where_angleInd_RS])**2)
                    PRS_f_fac = 1e23 * RForward2_RS * thicknessRS_forward / (Gamma[angleInd_RS+1] ** 2 * (1-beta[angleInd_RS+1]*cosAng[where_angleInd_RS])**2)

                    PRS_behind = PRS_b_fac * eats_function(pRS, nuPrimBehind, slow_coolingRS_behind, fast_coolingRS_behind, numRS[angleInd_RS], nucRS[angleInd_RS], PmaxF_RS[angleInd_RS], PmaxS_RS[angleInd_RS])
                    PRS_forward = PRS_f_fac * eats_function(pRS, nuPrimForward, slow_coolingRS_forward, fast_coolingRS_forward, numRS[angleInd_RS+1], nucRS[angleInd_RS+1], PmaxF_RS[angleInd_RS+1], PmaxS_RS[angleInd_RS+1])


                    PRSprimTemp = np.zeros(surfaceRings)
                    
                    PRSprimTemp[where_angleInd_RS] = (PRS_behind * weight1[where_angleInd_RS] + PRS_forward * weight2[where_angleInd_RS])  /  weight[where_angleInd_RS]
                    #PRSprimTemp =  PRS * 1e23 * R2 * thicknessRS / (Gamma[angleInd] ** 2 * (1-beta[angleInd]*cosAng)**2)

                if opticalDepth:
                    
                    alpha0Fforward = alpha0Ffactor * rhoPrimForward / B[angleInd+1]  *  gammac[angleInd+1]**(-5)
                    alpha0Fbehind = alpha0Ffactor * rhoPrimBehind / B[angleInd]  *  gammac[angleInd]**(-5)

                    alpha0Sforward =  alpha0Sfactor * rhoPrimForward  * gammam[angleInd+1]**(-5) /  B[angleInd+1]
                    alpha0Sbehind =  alpha0Sfactor * rhoPrimBehind  * gammam[angleInd]**(-5) /  B[angleInd]

                    numCentral = (num[angleInd] * weight1 + num[angleInd+1] * weight2) / weight
                    nucCentral = (nuc[angleInd] * weight1 + nuc[angleInd+1] * weight2) / weight
                    nuPrimCentral = (nuPrimBehind * weight1 + nuPrimForward * weight2) / weight
                    alpha0Fcentral = (alpha0Fbehind * weight1 + alpha0Fforward * weight2) / weight
                    alpha0Scentral = (alpha0Sbehind * weight1 + alpha0Sforward * weight2) / weight
                    alphanu = alpha0Fcentral * (numCentral>nucCentral) * ( (nuPrimCentral<=nucCentral) * (nuPrimCentral/nucCentral)**(-5/3.) + ((nucCentral<nuPrimCentral)&(nuPrimCentral<numCentral)) * (nuPrimCentral/nucCentral)**(-3) + (numCentral<=nuPrimCentral) * (numCentral/nucCentral)**-3*(nuPrimCentral/numCentral)**(-(p+5)/2.) )    +    alpha0Sforward * (numCentral<=nucCentral) * ( (nuPrimCentral<=numCentral)*(nuPrimCentral/numCentral)**(-5/3.) + ((numCentral<nuPrimCentral)&(nuPrimCentral<nucCentral)) * (nuPrimCentral/numCentral)**(-(p+4)/2.) + (nucCentral<=nuPrimCentral) * (nucCentral/numCentral)**(-(p+4)/2.) * (nuPrimCentral/nucCentral)**(-(p+5)/2.) )

                    tau = alphanu * (thickBehind * weight1 + thickForward * weight2) / 2 / weight     #tau is always negative
                    tau[tau<1e-8] = 0.
                    PprimTemp *= np.exp(tau)

                    
                    if reverseShock:
                        alpha0FRSforward = alpha0FRSfactor * rho3primForward / BRS[angleInd+1]  *  gammacRS[angleInd+1]**(-5)
                        alpha0FRSbehind = alpha0FRSfactor * rho3primBehind / BRS[angleInd]  *  gammacRS[angleInd]**(-5)

                        alpha0SRSforward =  alpha0SRSfactor * rho3primForward  * gammamRS[angleInd+1]**(-5) /  BRS[angleInd+1]
                        alpha0SRSbehind =  alpha0SRSfactor * rho3primBehind  * gammamRS[angleInd]**(-5) /  BRS[angleInd]

                        numRSCentral = (numRS[angleInd] * weight1 + numRS[angleInd+1] * weight2) / weight
                        nucRSCentral = (nucRS[angleInd] * weight1 + nucRS[angleInd+1] * weight2) / weight
                        alpha0FRScentral = (alpha0FRSbehind * weight1 + alpha0FRSforward * weight2) / weight
                        alpha0SRScentral = (alpha0SRSbehind * weight1 + alpha0SRSforward * weight2) / weight
                        alphanuRS = alpha0FRScentral * (numRSCentral>nucRSCentral) * ( (nuPrimCentral<=nucRSCentral) * (nuPrimCentral/nucRSCentral)**(-5/3.) + ((nucRSCentral<nuPrimCentral)&(nuPrimCentral<numRSCentral)) * (nuPrimCentral/nucRSCentral)**(-3) + (numRSCentral<=nuPrimCentral) * (numRSCentral/nucRSCentral)**-3*(nuPrimCentral/numRSCentral)**(-(pRS+5)/2.) )    +    alpha0SRSforward * (numRSCentral<=nucRSCentral) * ( (nuPrimCentral<=numRSCentral)*(nuPrimCentral/numRSCentral)**(-5/3.) + ((numRSCentral<nuPrimCentral)&(nuPrimCentral<nucRSCentral)) * (nuPrimCentral/numRSCentral)**(-(pRS+4)/2.) + (nucRSCentral<=nuPrimCentral) * (nucRSCentral/numRSCentral)**(-(pRS+4)/2.) * (nuPrimCentral/nucRSCentral)**(-(pRS+5)/2.) )
                        
                        tauRS = alphanuRS * (thickRSBehind * weight1 + thickRSForward * weight2) / 2 / weight     #tau is always negative
                        tauRS[tau<1e-8] = 0.
                        PRSprimTemp *= np.exp(tauRS + tau*2)


                
                angle_integ = (np.cos(xAngInt[:-1]) - np.cos(xAngInt[1:])) * phiInter   #Angular integration segments
                Pprim[rimI] = np.sum(PprimTemp  * angle_integ)     #In janskys
                if reverseShock: 
                    reverseComponent[rimI] = np.sum(PRSprimTemp * angle_integ)

                ### EATS integrator of thermal cocoon ###
                
                if thermalComp:# and (tobsRed[rimI] < 2e4):  #Equal arrival time surface integrator for thermal cocoon

                ### Finding distance corresponding to optical depth tau_photo
                    ### Need to find the elements closest to the EATSurface
                    near_ind_coc = np.argmin(np.abs((t_start + (R_photo_pre[:,0] - R0_coc) * tobs_conv_fac) * (1+z) - tobsRed[rimI])) ### The index in the precalculated arrays corresponding to the point closest time tobsRed[rimI] along the line of sight
                    near_ind_coc -= (t_start[near_ind_coc] + (R_photo_pre[near_ind_coc,0] - R0_coc) * tobs_conv_fac) * (1+z) > tobsRed[rimI]
                    cAngInt = np.linspace(0,upperLimitCoc[near_ind_coc],cocSurfaceRings+1) #Create surface rings
                    cAng = (cAngInt[:-1] + cAngInt[1:]) / 2     #Interpolate
                    cosAngCoc2 = np.cos(cAng)**2
                    nuPrimCoc = onePzFreq * GammaC * (1-beta_coc * np.cos(cAng))
                    use_thermal_here = np.ones(cocSurfaceRings,dtype=bool)
                    if near_ind_coc >= (precalc_length-1): 
                        use_thermal_here[0] = False
                        near_ind_coc -= 1
                    phiCInter = (cAng < lowerLimitCoc[near_ind_coc]) * 2*np.pi
                    if alpha != 0: #Off-axis
                        partialRingsIndC = np.where((cAng >= lowerLimitCoc[near_ind_coc]) & (cAng < upperLimitCoc[near_ind_coc]))   #Rings crossing the rim. Happens only when alpha != 0
                        offAxisFoVCoc = (theta_C[near_ind_coc]**2 - alpha**2 - cAng[partialRingsIndC]**2) / (2*alpha*cAng[partialRingsIndC])
                        offAxisFoVCoc[np.where(offAxisFoVCoc<-1)] = -1.
                        if np.sum(offAxisFoVCoc > 1) != 0: offAxisFoVCoc[np.where(offAxisFoVCoc > 1)] = 1.
                        phiCInter[partialRingsIndC] = 2*np.pi - 2*np.arccos(offAxisFoVCoc)


                    for i_coc_rings in range(cocSurfaceRings):
                        if i_coc_rings != 0:
                            ph_tobs_below = False
                            sign_change = 0
                            while True:
                                tobs_c_below = (t_start[near_ind_coc] + (R_photo_pre[near_ind_coc,i_coc_rings] - R0_coc) * tobs_conv_fac) * (1+z) 
                                if tobs_c_below > tobsRed[rimI]: 
                                    near_ind_coc -= 1
                                    if sign_change == -1: break
                                    else:
                                        if near_ind_coc <= 0: break
                                        sign_change = 1
                                else:
                                    near_ind_coc += 1
                                    if near_ind_coc >= (precalc_length-2): 
                                        near_ind_coc -= 1
                                        use_thermal_here[i_coc_rings] = False
                                        break
                                    if sign_change == 1: break
                        else: tobs_c_below = (t_start[near_ind_coc] + (R_photo_pre[near_ind_coc,0] - R0_coc) * tobs_conv_fac) * (1+z)
                        if use_thermal_here[i_coc_rings]:

                            tobs_c_above = (t_start[near_ind_coc+1] + (R_photo_pre[near_ind_coc+1,i_coc_rings] - R0_coc) * tobs_conv_fac) * (1+z)
                            weight_c1 = tobsRed[rimI] - tobs_c_below
                            weight_c2 = tobs_c_above - tobsRed[rimI]
                            weight_c = tobs_c_above - tobs_c_below
                            
                                        
                        ### Interpolating

                            R_photo_int[i_coc_rings] = (R_photo_pre[near_ind_coc,i_coc_rings] * weight_c2 + R_photo_pre[near_ind_coc+1,i_coc_rings] * weight_c1) / weight_c
                            R_coc_ind = np.argmin(np.abs(R_coc - R_photo_int[i_coc_rings]))  ### What index in the precalculated arrays corresponds to the photospheric radius?
                            R_coc_ind -= R_coc[R_coc_ind] < R_photo_int[i_coc_rings]
                            T_int[i_coc_rings] = (T[R_coc_ind] * (R_coc[R_coc_ind+1]-R_photo_int[i_coc_rings]) + T[R_coc_ind+1] * (R_photo_int[i_coc_rings]-R_coc[R_coc_ind])) / (R_coc[R_coc_ind+1] - R_coc[R_coc_ind])

              

                    thermal_elements = np.where(use_thermal_here)                    
                    hkBT_coc = hcgs / kB / T_int[thermal_elements]

                    thermalExpTerm = nuPrimCoc[thermal_elements]*(1+z)*hkBT_coc

                    if np.sum(thermal_elements)==0: thermal_component[rimI] = 0.
                    else:
                        if (np.min(thermalExpTerm) < -200) or (np.max(thermalExpTerm) > 200): thermal_component[rimI] = 0. 
                        else:
#                        thermalForward = 2 * hcgs * (nuPrimForwardCoc*(1+z))**3 / (np.exp(thermalExpTermF) - 1) * np.cos(cAng)
#                        thermalBehind = 2 * hcgs * (nuPrimBehindCoc*(1+z))**3 / (np.exp(thermalExpTermB) - 1) * np.cos(cAng)
                       
#                        thermalForward  *= Rcoc[angCocInd+1] ** 2 / (GammaC ** 2 * (1-beta_coc[angCocInd+1]*cosAngCoc2)**2)
#                        thermalBehind *= Rcoc[angCocInd] ** 2 / (GammaC ** 2 * (1-beta_coc[angCocInd]*cosAngCoc2)**2)
                            thermal_component[rimI] =  np.sum( 2 * hcgs * (nuPrimCoc*(1+z))**3 / (np.exp(thermalExpTerm) - 1) * np.cos(cAng[use_thermal_here]) * R_photo_int[use_thermal_here] ** 2 / (GammaC ** 2 * (1-beta_coc*cosAngCoc2[use_thermal_here])**2))
            






         
                    
                    #                rho_coc[0] = np.copy(rho_coc[1])
                
                #plt.show()

                #Setting these two first lines to see evaluation


                
                #plt.show()


#                    indRimCoc = np.argmin(np.abs(tobsRimNow-tobsRimC))    #Find what index the point at the rim has with the observer's time we are looking for.
#                    indRimCoc -= (tobsRimC[indRimCoc] > tobsRimNow)  

#                    if indRimCoc >= kill_cocoon:   #If the optical depth has decreases to make the cocoon transparent and the observer time of the rim has passed this
#                        thermal_component[rimI] = 0.
#                    else:    

                        #Evaluating the field of view. The field of view has crossed the rim of the cocoon if an element of fovInRimCoc == False.
                    #upLim1 , upLim2 = (fovInRimCoc[indRimCoc]==False) , (theta_C[indRimCoc]+alpha) , 
                    #upperLimitCoc = upLim1 * upLim2 + fovInRimCoc[indRimCoc]*fovAngleCoc[indRimCoc]
                    #lowerLimitCoc = (fovInRimCoc[indRimCoc]==False)*(theta_C[indRimCoc]-alpha) + fovInRimCoc[indRimCoc]*fovAngleCoc[indRimCoc]

                    #cAngInt = np.linspace(0,upperLimit,cocSurfaceRings+1) #Create surface rings
                    #cAng = (cAngInt[:-1] + cAngInt[1:]) / 2     #Interpolate


                    ### Finding what ring segment to use - if photosphere==True the angle between the observer and the normal of the ring segment needs to be taken into account by correcting after the normal index finding process
                    
                    #angCocInd[0] = np.argmin(np.abs(tobsC-tobsRed[rimI]))
                        
                    #    for cocAng in range(1,cocSurfaceRings): #Finding index to each ring segment
                    #        angCocInd[cocAng] = np.copy(angCocInd[cocAng-1])
                    #        while tobsRed[rimI] < ((1+z) * (tburstC[angCocInd[cocAng]] - (Rcoc[angCocInd[cocAng]]-R0_C) * np.cos(cAng[cocAng])/c)):
                    #            angCocInd[cocAng] -= 1
                    #        if photosphere: #Correcting for photosphere. See note from 2015-03-05
                            
                            
                    #            Rcoc_photo_array = Rcoc - Dtau * np.cos(cAngInt[cocAng])   ### Which Rcoc corresponds to the photospheric radius equaling the radius that the EATS integrator wants. See note from 2015-03-05
                            
                    #            Rcoc_corr_ind = np.argmin(np.abs(Rcoc_photo_array - Rcoc[angCocInd[cocAng]]))   #Assigning index corresponding to the radial step with a photosphere corresponding to the radius that the EATS integrator is looking for
                    #            Rcoc_corr_ind -= ((Rcoc[Rcoc_corr_ind] - Dtau[Rcoc_corr_ind] * np.cos(cAngInt[cocAng])) > Rcoc[angCocInd[cocAng]])  #Correcting to make sure we are behind the sought index
                    #            angCocInd[cocAng] = np.copy(Rcoc_corr_ind)
                        
                    #not_transparent = angCocInd < kill_cocoon




                   #     tobsCBehind = (1+z) * (tburstC[angCocInd] - (Rcoc[angCocInd]-R0_C) * np.cos(cAng) / c)
                   #     try:tobsCFront = (1+z) * (tburstC[angCocInd+1] - (Rcoc[angCocInd+1]-R0_C) * np.cos(cAng) / c)
                   #     except:
                   #         print angCocInd
                   #         print len(tburstC)
                   #         print 'Crashed!'
                    
                   #     weightC,weightC1,weightC2 = np.log10( tobsCFront / tobsCBehind )  ,  np.log10( tobsCFront / tobsRimNow )  ,  np.log10( tobsRimNow / tobsCBehind )
                            
                    
            
        

        
        if (not RS_only) and (not thermal_only):
            F = (1+z)/(2*D**2)*Pprim


        if reverseShock: 
            if RS_only: F = (1+z)/(2*D**2)*reverseComponent #If RS_only is true, only reverse shock component will be plotted
            elif not thermal_only: 
                FRS = (1+z)/(2*D**2)*reverseComponent
                
                if plotComponents & (runOption == 'LC') & (not FS_only):    #Plot RS and FS components separately?
                    colourCycle = ['b','g','r','c','m','y','k']     #Cycle of colours in plotting. Matplotlib standard cycle
                    scalePlotTime = {'d': 86400. , 'h' : 3600. , 'm' : 60. , 's' : 1. }
                    scaleFluxAxis = {'mJy' : 1.e3 , 'Jy' : 1. }
                    
                    plt.plot(tobsRed/scalePlotTime[daysOrSec],F * scaleFluxAxis[fluxAxis],'%s--'%colourCycle[nuIte%len(colourCycle)])
                    plt.plot(tobsRed/scalePlotTime[daysOrSec],FRS * scaleFluxAxis[fluxAxis],'%s-.'%colourCycle[nuIte%len(colourCycle)])
                F += FRS    #If FS_only is true, only forwards shock component will be plotted

        if thermalComp:
            Fthermal = (1+z)/(2*D**2)*thermal_component
            if plotComponents and (runOption == 'LC'):
                colourCycle = ['b','g','r','c','m','y','k']     #Cycle of colours in plotting. Matplotlib standard cycle                                                
                scalePlotTime = {'d': 86400. , 'h' : 3600. , 'm' : 60. , 's' : 1. }
                scaleFluxAxis = {'mJy' : 1.e3 , 'Jy' : 1. }
                if not reverseShock and not thermal_only: plt.plot(tobsRed/scalePlotTime[daysOrSec],F * scaleFluxAxis[fluxAxis],'%s--'%colourCycle[nuIte%len(colourCycle)])
                plt.plot(tobsRed/scalePlotTime[daysOrSec],Fthermal * scaleFluxAxis[fluxAxis],'%s:'%colourCycle[nuIte%len(colourCycle)])
            if (not RS_only) and (not FS_only):
                if thermal_only: F = Fthermal
                else:  F += Fthermal

        if noChi2 == False: 
            if runOption=='LC':
                for fInd2 in range(numberOfEmpties[nuIte]):
                    #testF[fInd2] = F[np.argmin(np.abs(tdataLC[fInd2]-tobsRed))]
                    middleInd = np.argmin(np.abs(tdataLC[fInd2]-tobsRed))
                    behindInd = middleInd - (tdataLC[fInd2] < tobsRed[middleInd])
 
                    Fweight1 , Fweight2 , Fweight = np.log10(tdataLC[fInd2]) - np.log10(tobsRed[behindInd]) , np.log10(tobsRed[behindInd+1]) - np.log10(tdataLC[fInd2]) , np.log10(tobsRed[behindInd+1]) - np.log10(tobsRed[behindInd])
                    Finter = (F[behindInd] * Fweight2 + F[behindInd + 1] * Fweight1) / Fweight

                    if chi2_type == 'lin': chi2 += ((Fdata[fInd2] - Finter) / errorbar[fInd2])**2
                    elif chi2_type == 'log': chi2 += (np.log10(Fdata[fInd2]/Finter)) ** 2 / np.log10((Fdata[fInd2]+errorbar[fInd2])/Fdata[fInd2])**2
                    else: 
                        print 'Bad chi2 type %s. Now exiting'%chi2_type
                        raise SystemExit(0)
                
            else:
                #if (np.sum(F < 0) == 0): # Avoiding negativ fluxes
                
                if chi2_type == 'lin': 
                    F[np.where(F<0)] = 0.
                    chi2 += np.sum(((Fdata - F) / errorbar)**2)
                elif chi2_type == 'log': 
                    F[np.where(F<=0)] = 1e-30
                    chi2 += np.sum(np.log10(Fdata/F)**2 / np.log10((Fdata+errorbar)/Fdata)**2)
                else: 
                    print 'Bad chi2 type %s. Now exiting'%chi2_type
                    raise SystemExit(0)
                #else:
                #    print 'Bad flux output. Returning chi2 = \'inf\''
                #    return float('inf'),float('nan'),float('nan') #If we get a negativ flux, return 'NaN'
        
        if (runOption == 'LC') or (runOption == 'one-sigma'): 
#            print "Loop time = %f"%loopTimerTotal
            lightcurve[nuIte] = np.copy(F)#np.concatenate([F ,   [-1]*(len(lightcurve[nuIte]) - len(F))])

    
        
    if runOption == 'LC': 
        print "Synchrotron time use: %f s"%(time.time() - dynTimeStart)

        if allowPrint & (noChi2 == False): print "chi2 = %s\nReduced chi2 = %s"%(chi2,chi2 / (numberOfPoints-ndims) )
        #print lightcurve

        
                #outputCounter = 0
            #else: outputCounter += 1

        startJB_index = np.argmin(np.abs(theta-alpha-fovAngle))

        if startJB_index >= len(tobs)-1: startJetBreak = - tobs[startJB_index] ### If jetbreak starts at really late times
        else:
            startJB_index -= ((theta[startJB_index] - alpha - fovAngle[startJB_index]) < 0)  ### Making sure the field of view is still just a little bit smaller than the rim
            startJB_weight1 = theta[startJB_index] - alpha - fovAngle[startJB_index]
            startJB_weight2 = fovAngle[startJB_index+1] - theta[startJB_index+1] + alpha
            print startJB_weight1
            print startJB_weight2
        
        endJB_index = np.argmin(np.abs(theta+alpha-fovAngle))
        if endJB_index >= len(tobs)-1: endJetBreak = - tobs[endJB_index]  ### if jetbreak starts at really late times
        else:    

            endJB_index -= ((theta[endJB_index] + alpha) < fovAngle[endJB_index])  ### Making sure the field of view is still just a little bit smaller than the last crossing of the rim
            endJB_weight1 = theta[endJB_index] + alpha - fovAngle[endJB_index]
            endJB_weight2 = fovAngle[endJB_index+1] - theta[endJB_index+1] - alpha


        if startJB_index < len(tobs)-1: startJetBreak = (tobs[startJB_index] * startJB_weight2 + tobs[startJB_index+1] * startJB_weight1) / (startJB_weight1 + startJB_weight2)
        if endJB_index < len(tobs)-1: endJetBreak = (tobs[endJB_index] * endJB_weight2 + tobs[endJB_index+1] * endJB_weight1) / (endJB_weight1 + endJB_weight2)

        startJB_text = '%s %f'%('='*(startJetBreak > 0) + '>'*(startJetBreak < 0) , startJetBreak*(((startJetBreak>0)*2) - 1) / 86400)
        endJB_text = '%s %f'%('='*(endJetBreak > 0) + '>'*(endJetBreak < 0) , endJetBreak*(((endJetBreak>0)*2) - 1) / 86400)

        print "Field of view started crossing the rim at tobs %s days and covered the entire rim at tobs %s days."%(startJB_text,endJB_text)

        return lightcurve , startJetBreak , endJetBreak , tobsGrid
    elif runOption == 'one-sigma': return lightcurve , 0. , 0.  ,  tobsGrid
    else: 
        
        return chi2 , None , None
