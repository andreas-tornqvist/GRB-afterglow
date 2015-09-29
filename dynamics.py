#This module considers the Nava et al. (2013) paper. 
def dyn(R,Gamma0,epsilone,epsilonp,epsilone3,epsilonp3,eB,eB3,E0,M0,M03,s,A0,t0,p,pRS,theta0,z,tprompt,reverseShock,tobsEnd,profile_cutoff):
    import numpy as np
    import matplotlib.pyplot as plt
    c,pi,mp,me,kB = 2.9979e10,np.pi, 1.6726e-24, 9.1094e-28 , 1.38065e-16
    qe = 4.803204e-10
    sigma_B = 5.6704e-5  ### Stephan Boltzann constant in cgs                   
    rad_const = 4 * sigma_B / c   #### Radiation constant                 
#Initial conditions
    #M0 = E0*c**-2
    if reverseShock: grandNan = [float('nan')] * 14
    else: grandNan = [float('nan')] * 9
    startSpread = 1e12
    useSpread = True   #If you want to turn spread off, turn it off on this line! (or the next)
    printStats = False
    
    RScutOff = 0.
    
    mup,mue = 1.,me / mp
    mmean = (me+mp)/2
    sigmaT = 6.6524e-25
    
    
    zero = np.concatenate([np.array([1.]),np.zeros(len(R)-1)])
    allZero = np.zeros(len(R))
    beta0 = (1 - 1/Gamma0**2)**(1/2.)
    beta = beta0*zero
    V2 = np.zeros(len(R))

    
    
    if A0 == 0: #Cocoon
        #T0 = epsilone
        tburst , tobs , tcomoving = np.concatenate([np.array([1.]),np.zeros(len(R)-1)]) ,np.concatenate([np.array([1.]),np.zeros(len(R)-1)])  , np.concatenate([np.array([1.]),np.zeros(len(R)-1)]) 
        rho = np.zeros(len(R))
        nprim = np.zeros(len(R))
        rho04 = np.zeros(len(R))
    else:
        tburst,tobs,tcomoving = zero* R[0]/c/beta0 , zero*R[0]/c*(1./beta0 - 1)*(1+z) , zero *  R[0]/c/beta0/Gamma0     #zero* R[0]/(1.-1./Gamma[0]**2)**(1/2.)/c ,
        rho = A0 * R**(-s) #* CBM76 * 0. #OBS! This may only be used for the energy evaluation!
        rhoprim = np.zeros(len(R))
        if (s != 0) and (R[-1] > 10**profile_cutoff): rho[np.where(R>10**profile_cutoff)] = rho[np.argmin(np.abs(R-10**profile_cutoff))]
        rho04 = 3 * M0 / (2*pi*(1-np.cos(theta0)) * (  c**3*beta0**3 * (tburst[0]**3 - tprompt**3)  ) )
     

    theta = zero * theta0 #if useSpread == False, this value is fixed, otherwise it changes
#    theta[0] = np.pi/2
    Delta0 = tprompt
    m0 = 2*pi*R[0]**3*(1-np.cos(theta[0])) * rho[0] /3

    eiFrac = .8 #Fraction of the swept-up mass left to sweep up. Determines the initial fraction of mass in shock front
    checkLength = 1-eiFrac   #Test to see if energy injection from prompt outflow integration works
    
    
        #Setting initial values
    gammaAdi,dGammaEffdGamma, = zero*(4+1/Gamma0)/3.,zero*0. #Good guess for dGammaEffdGamma? (coasting phase)
    Gamma,GammaEff,m = zero*Gamma0,zero*(gammaAdi[0]*Gamma0**2-gammaAdi[0]+1)/Gamma0,zero*m0
    Bconst = 32.*pi*mp
    Eint2 =   zero*0.
    Bterm_array = np.copy(zero)

    EForPlot = zero*0.
    epsilonSave = zero
    if reverseShock:
        rho4 =  zero*M0 / (4*pi*R[0]**2*Delta0) 
        gamma43 = zero *  Gamma0**2*(1 - beta[0]**2)
        if gamma43[0] <= 1.: gamma43[0] = 1.+1e-15
        n3prim = rho4[0] * 4 * gamma43[0]
        Eint3 = zero*0.
        gammaAdi3,dGammaEff3dGamma = zero*(4+1/Gamma0)/3 , zero*0
        GammaEff3 = zero*(gammaAdi3[0]*Gamma0**2-gammaAdi3[0]+1)/Gamma0
        beta43 = np.copy(zero)
        V3 = np.zeros(len(R))
    if A0==0: 
        T = zero
        R_photo = zero
        rho_coc = np.zeros(len(R))
        Rcoc_dist = np.zeros(len(R))
        Rcoc_dist[0] = (R[1] - R[0]) / 2
        Rcoc_dist[-1] = (R[-1] - R[-2]) / 2
        Rcoc_dist[1:-1] = (R[2:] - R[:-2]) / 2
        tau_coc = 1  ### At which optical depth the photosphere exists 

    M3  =   zero*(1-eiFrac) * M0            #M0 *(beta0-beta[0])/beta[0]   /   Delta0   * R[0]
    dM3dR = 0.
    shutOff = False             #If the reverse shock feed has shut off or not. Used for printing of statistics only
    cutOffSlope = 0.1
    
    i = 0
    while i < len(R):#for i in xrange(len(R)-1):
        if A0==0:print i
        #dmdR = 2*pi*R[i]**(2.-s)*A0*(1-np.cos(theta[i])) #Will differ with different geometry. Make sure that theta is defined here
        #dmdR = 
            #m[i+1] = 4*pi/(3-s)*R[i]**(3-s)*A0#
        #m[i+1] = m[i] + dmdR*(R[i+1]-R[i]) #Numeric stepper. Defined value will not be used until next iteration! m[0]=0
        dmdR = 2*pi*R[i]**(2.-s)*A0*(1-np.cos(theta[i])) #Will differ with different geometry. Make sure that theta is defined here
	if s==0: 
            
            m[i] = rho[i] * 2*pi*R[i]**3/3*(1-np.cos(theta[i]))
	else: 
            delta_m = 2*pi*(R[i+1]**3-R[i]**3)*(rho[i+1]+rho[i])*(1-np.cos(theta[i]))/6 #Factor 6 is from averaging rho and dividing by 3 for the volume
            m[i+1] = m[i] + delta_m
            #dmdR = delta_m / (R[i+1] - R[i])
        
            #m[i+1] = m[i] + dmdR*(R[i+1]-R[i]) 

        if useSpread:
            vperp = ( gammaAdi[i]*(gammaAdi[i]-1)*(Gamma[i]-1) / (1 + gammaAdi[i]*(Gamma[i]-1)) )**(1/2.) * c
        if i == 0: 
            Rdist = np.array([(R[1]-R[0])/2.])
            dthetadR = np.array([0.])
#            V2[0] = m[0] / rho[0] / 4/ Gamma0
                #thetaDist = np.array([(theta[1]-theta[0])/2.]) #Until we assume an angular distribution, we can integrate over opening angle analytically
        else: 
            Rdist = np.append(Rdist,(R[i]-R[i-1])/2.) #The distance between each radial step. The second half added in the end of the loop. For the integrator
            if useSpread:
                dthetadR = np.append(dthetadR,2*vperp/(R[i]+R[i+1])/Gamma[i]/beta[i]/c)               ###(theta[i]-theta[i-1])/(R[i]-R[i-1]))
            else:
                dthetadR = np.append(dthetadR,0.)
                
                #thetaDist = np.append(thetaDist,(theta[i]-theta[i-1])/2.) #Make sure that theta[i] is defined.
            #pp2,pe2 = ((Gamma[:i+1]-1)*epsilonp/mup + 1)**2 -1 , ((Gamma[:i+1]-1)*epsilone/mue + 1)**2 - 1

        

        #Radiative efficiency estimation

        if i == 0: 
            epsilon = np.copy(zero)*0.
            epsilonRS = np.copy(zero)*0.
            if A0 == 0: #Cocoon
                r0 = 1e7
                T[0] = (E0 / 4 / np.pi / r0**2 / c / rad_const) ** 0.25 / Gamma[0] * (R[0] / r0) ** (-2/3.)
                gammae = zero * T[0] * kB * 3 / 2 *c **-2 / me + 1
                #print gammae[0]
                gammap = mue / mup * (gammae-1) + 1
#                gammap,gammae = zero * (Gamma[i]-1)*mue + 1 , zero * (Gamma[i]-1)*mup + 1
                pp2,pe2 = zero * gammap[0]**2 -1 , zero * gammae[0]**2 - 1 #Interesting thought: this equation is the only place where the radiative efficiency epsilon actually does matter.
                Eint2[0] = M0 * (mup * (gammap[0] - 1) + mue * (gammae[0] - 1)) * c**2
                print 'Eint2[0] = %s'%Eint2[0]
#                T = zero * T0#Eint2[0] * (me+mp) / M0 / 3 / kB
            else:
                gammap,gammae = (Gamma-1)*epsilonp/mup + 1 , (1-epsilon)*(Gamma-1)*epsilone/mue + 1
                pp2,pe2 = gammap**2 -1 , gammae**2 - 1 
                if reverseShock:
                    gammapRS , gammaeRS = 1 + epsilonp3 * (gamma43 - 1) , 1 + epsilone3*(gamma43 - 1) * mp/me
                    pp2rs,pe2rs = gammapRS**2 - 1   ,   gammaeRS**2 - 1
        if True:
            if (A0 != 0):    # Jet dynamics
                gammap[i],gammae[i] = (Gamma[i]-1)*epsilonp/mup + 1 , (1-epsilon[i])*(Gamma[i]-1)*epsilone/mue + 1
                nprim = rho[i] / mp * 4 * Gamma[i]    #Shock front density 
                
                    
                B = (Bconst*eB*Gamma[i]*(gammae[i]-1)*nprim)**(1/2.)*c
                
                ### Radiated energy estimation. See notes from 8/9 -15
                pimec = 2 * pi * me * c
                gamToNuFactor = qe * B / pimec

                gammamFac1 = (p-2)/(p-1)
                gammamFac2 = epsilone*mp/me
                gammam = gammamFac1 * (gammamFac2 * (GammaEff[i] -1) + 1)
                gammac1t = 3*pimec/sigmaT/tobs[i]
                gammac = gammac1t/Gamma[i]/B**2

                num = gammam**2 * gamToNuFactor
                nuc = gammac**2 * gamToNuFactor
                
                phipF = 0.54 + 0.08*p
                phipS = 1.89 - 0.935*p + 0.017*p**2

                PmaxF = phipF * 2.234*qe**3*B/me/c**2
                PmaxS = phipS*11.17*(p-1)*qe**3*B/(3*p-1)/me/c**2

                if num < nuc:   #Slow cooling regime
#                    P_total = PmaxS * (3/4. + 2/(p-3) * (1-(nuc/num)**(-(p-3)/2)) + 2/(p-2) * (nuc/num)**(-(p-1)/2))
                    P_total_term = (3/4. + 2/(p-3) * (1-(nuc/num)**(-(p-3)/2)) + 2/(p-2) * (nuc/num)**(-(p-1)/2))
                    P_total = PmaxS * P_total_term
                    
                else:  # Fast cooling regime
                    P_total = PmaxF * (3/4. + 2*((num/nuc)**(0.5) - 1) + (num/nuc)**-.5 * (2/(p-2)))
                P_total *= m[i] / 2 / (me+mp)
                dEsh_dt_co = Gamma[i]*(Gamma[i]-1)*beta[i] * dmdR * c**3
                
                if dEsh_dt_co == 0.: epsilon[i+1] = 0.
                else: epsilon[i+1] = P_total / dEsh_dt_co

                if reverseShock:
                    n3prim = rho4[i] * 4 * gamma43[i]
                    gammapRS[i] , gammaeRS[i] = 1 + epsilonp3 * (gamma43[i] - 1) , 1 + (1-epsilonRS[i])*epsilone3*(gamma43[i] - 1) * mp/me
                    pp2rs[i],pe2rs[i] = gammapRS[i]**2 - 1   ,   gammaeRS[i]**2 - 1
                    BRS = (Bconst*eB3*n3prim*Gamma[i]*(gammaeRS[i]-1))**(1/2.)*c
                    if BRS != 0:
                        gamToNuFactorRS = qe * BRS / pimec

                        gammamFac1RS = (pRS-2)/(pRS-1)
                        gammamFac2RS = epsilone3*mp/me
                        gammamRS = gammamFac1RS * (gammamFac2RS * (GammaEff3[i] -1) + 1)
                        gammac1tRS = 3*pimec/sigmaT/tobs[i]
                        gammacRS = gammac1tRS/gamma43[i]/BRS**2

                        numRS = gammamRS**2 * gamToNuFactorRS
                        nucRS = gammacRS**2 * gamToNuFactorRS
                
                        phipFRS = 0.54 + 0.08*pRS
                        phipSRS = 1.89 - 0.935*pRS + 0.017*pRS**2

                        PmaxFRS = phipFRS * 2.234*qe**3*BRS/me/c**2
                        PmaxSRS = phipSRS*11.17*(pRS-1)*qe**3*BRS/(3*pRS-1)/me/c**2

                        if numRS < nucRS:   #Slow cooling regime

                            P_total_termRS = (3/4. + 2/(pRS-3) * (1-(nucRS/numRS)**(-(pRS-3)/2)) + 2/(pRS-2) * (nucRS/numRS)**(-(pRS-1)/2))
                            P_totalRS = PmaxSRS * P_total_termRS
                    
                        else:  # Fast cooling regime
                            P_totalRS = PmaxFRS * (3/4. + 2*((numRS/nucRS)**(0.5) - 1) + (numRS/nucRS)**-.5 * (2/(pRS-2)))
                        P_totalRS *= M3[i] / 2 / (me+mp)
                        dEsh_dt_coRS = Gamma[i]*(gamma43[i]-1)*beta[i] * dM3dR * c**3
                        if dEsh_dt_coRS == 0: epsilonRS[i+1] = 0.
                        else:

                
                            epsilonRS[i+1] = P_totalRS / dEsh_dt_coRS
                            if np.abs(epsilonRS[i+1]) < 1e-10: epsilonRS[i+1] = 0
                            if epsilonRS[i+1] < 0:
                                epsilonRS[i+1] = 0.
                    else:
                        epsilonRS[i+1] = 0.
                if epsilonRS[i+1] >= 1.: epsilonRS[i+1] = 1.-1e-5
            if epsilon[i+1] >= 1.: epsilon[i+1] = 1.-1e-5


            pp2[i],pe2[i] = gammap[i]**2 -1 , gammae[i]**2 - 1 #Interesting thought: this equation is the only place where the radiative efficiency epsilon actually does matter.
            
       #We do not need to include the thickness - we know how many particles are radiating at every moment!
            
            
                
        
            
            
        if A0 == 0:  #Calculating cocoon
            
            #Eint2[i] = M0 * (mup * (gammap - 1) + mue * (gammae - 1))
            #print "gammap = %s"%gammap[i]
            #print pp2[i]
            #print pe2[i]
            #raw_input(gammae[0])
            Bterm = M0 * (mup * pp2[0] / gammap[0] + mue * pe2[0] / gammae[0])
            dGammadR = GammaEff[i]/R[i] * Bterm  /  (M0 + Eint2[i] * dGammaEffdGamma[i] - GammaEff[i] / 3 / Gamma[i] * Bterm)
            #print dGammadR
            Gamma[i+1] = Gamma[i] + dGammadR * (R[i+1] - R[i])
            print Gamma[i+1]
        else: #Calculating jet
            
            Rfraq = (R[:i+1]/R[i])**2
            Gammafraq = (Gamma[i]/Gamma[:i+1])**(2/3.)
            RGammaFraqP = Rfraq*Gammafraq*pp2[:i+1]+1
            RGammaFraqE = Rfraq*Gammafraq*pe2[:i+1]+1

            RGFraqSqrtP = (RGammaFraqP)**(1/2.)
            RGFraqSqrtE = (RGammaFraqE)**(1/2.)
            
            RdistR2 = Rdist * R[:i+1]**2
            Rsum2 = R[:i+1]**2
            ### New adiabatic expansion approach - see note from 27/8-15
            
            Eint2[i] = -2*pi*c**2*np.sum(RdistR2*rho[:i+1]*(mup*(RGFraqSqrtP-1) + mue*(RGFraqSqrtE-1)) * (np.cos(theta[:i+1])-1)) #Array multiplicative integrator
            if R[i] > profile_cutoff: dlnrhodR = 0.
            else: dlnrhodR = s/R[i]
            rhoprim[i] = 4*Gamma[i]*rho[i]
            V2[i] = m[i] / rhoprim[i]

            f_2 = GammaEff[i]*(gammaAdi[i]-1)*Eint2[i]/V2[i]*m[i]/Gamma[i]/rhoprim[i]

            
            
                


            if reverseShock:
                    
#                beta[i] = (1 - 1/Gamma[i]**2) ** (1/2.)

                    
                #Pressure terms
                

                    #Approach: The region 2 (Nava 2013) is the forward shock, i.e. same as in the case without the reverse shock. Now we are adding terms for the reverse shocks, in the region labeled 3.
                
                
                RGammaFraqPrs = Rfraq*Gammafraq*pp2rs[:i+1]+1
                RGammaFraqErs = Rfraq*Gammafraq*pe2rs[:i+1]+1

                RGFraqSqrtPrs = (RGammaFraqPrs)**(1/2.)
                RGFraqSqrtErs = (RGammaFraqErs)**(1/2.)

                Eint3[i] = -2*pi*c**2*np.sum(RdistR2*rho[:i+1]*(mup*(RGFraqSqrtPrs-1) + mue*(RGFraqSqrtErs-1)) * (np.cos(theta[:i+1])-1)) #Array multiplicative integrator

                if n3prim == 0:
                    f_3 = 0.
                else:
                    V3[i] = M3[i] / n3prim / mmean
                    f_3 = GammaEff3[i]*(gammaAdi3[i]-1)*Eint3[i]/V3[i]*M3[i]/gamma43[i]/rhoprim[i] * (0.5/Gamma0 - Gamma0 / 2/Gamma[i]**2)
                
#                Eint2[i] = 2*pi*c**2*np.sum(RdistR2* ( dGammaEffdGamma[i] *  rho[:i+1]*(mup*(RGFraqSqrtP-1) + mue*(RGFraqSqrtE-1))   +   dGammaEff3dGamma[i] *  rho4[i] * (beta0-beta[:i+1])/beta[:i+1] * (mup*((RGammaFraqPrs)**(1/2.)-1) + mue*((RGammaFraqErs)**(1/2.)-1))      ) * (1 - np.cos(theta[:i+1]))) #Array multiplicative integrator
                
#                dEaddRterm = -np.sum(Rsum2*  ( GammaEff[i] *  rho[:i+1]*Rfraq * Gammafraq*(mup*pp2 / RGFraqSqrtP + mue*pe2/RGFraqSqrtE)   +     GammaEff3[i] *  rho4[i] * (beta0 - beta[:i+1])/beta[:i+1]   *     Rfraq * Gammafraq*(mup*pp2rs / (RGammaFraqPrs)**(1/2.) + mue*pe2rs/(RGammaFraqErs)**(1/2.))        ) * (np.cos(theta[:i+1])-1) * Rdist) #Note from 16/8 -13 with RS addition from 28-29/1 -14
                
#                Bterm = np.sum((GammaEff[i] * rho[:i+1]*(mup*RGFraqSqrtP + mue*RGFraqSqrtE)  +  (beta0-beta[:i+1])/beta[:i+1]   *   GammaEff3[i]  *      rho4[i] * (mup*(RGammaFraqPrs)**(1/2.) + mue*(RGammaFraqErs)**(1/2.))   )*np.sin(theta[:i+1])*dthetadR*Rsum2*Rdist)
#                print np.shape(dmdR)
#                print np.shape(dM3dR)
#                print np.shape(f_2)
#                print np.shape(f_3)
#                print np.shape(Gamma[i])
#                print np.shape(Eint3[i])
#                print np.shape(dlnrhodR)
#                print np.shape(m[i])
#                print np.shape(gammaAdi3[i])
#                print np.shape(M3)
#                print np.shape(((Gamma[i]-1)*(GammaEff[i]+1)*dmdR*c**2 - GammaEff[i]*(gammaAdi[i]-1)*Eint2[i]*(dmdR/m[i] + dlnrhodR) - GammaEff3[i]*(gammaAdi3[i]-1)*Eint3[i]*(dM3dR/M3[i] - 2/R[i])) / ((M3+m[i])*c**2 + Eint2[i]*dGammaEffdGamma[i] + Eint3[i]*dGammaEff3dGamma[i] + f_2 + f_3) * (R[i+1] - R[i]))
                
                Gamma[i+1] = Gamma[i] - ((Gamma[i]-1)*(GammaEff[i]+1)*dmdR*c**2 - GammaEff[i]*(gammaAdi[i]-1)*Eint2[i]*(dmdR/m[i] + dlnrhodR) - GammaEff3[i]*(gammaAdi3[i]-1)*Eint3[i]*(dM3dR/M3[i] - 2/R[i])) / ((M3[i]+m[i])*c**2 + Eint2[i]*dGammaEffdGamma[i] + Eint3[i]*dGammaEff3dGamma[i] + f_2 + f_3) * (R[i+1] - R[i])
                if np.isnan(Gamma[i+1]):
                    print 'isnan!'
                if Gamma[i+1] < 1:
                    Gamma[i+1] = 1.+1e-4
#                    plt.plot(R[:i+1], Gamma[:i+1])
#                    plt.loglog()
#                    plt.show()

                ##### TO DO ######
            

            
#                Gamma[i+1] = Gamma[i] - ((GammaEff[i]+1)*(Gamma[i]-1)*c**2*dmdR + (Gamma[i]-Gamma0+GammaEff3[i]*(gamma43[i]-1)) * dM3dR * c**2         - (2*pi*c**2/R[i] + np.sin(theta[i])/(1-np.cos(theta[i]))*dthetadR[i] )*dEaddRterm - 2*pi*c**2*Bterm) / ((M3[i]+m[i])*c**2 + Eint2[i] + 4*pi*c**2/3./Gamma[i] * dEaddRterm) * (R[i+1] - R[i]) #Note from 25/8-13
            else:  #Without reverse shock
                
                Gamma[i+1] = Gamma[i] - ((Gamma[i]-1)*(GammaEff[i]+1)*dmdR*c**2 - GammaEff[i]*(gammaAdi[i]-1)*Eint2[i]*(dmdR/m[i] + dlnrhodR)) / ((M0+m[i])*c**2 + Eint2[i]*dGammaEffdGamma[i] + f_2) * (R[i+1] - R[i])
                                         
                if ((Gamma[i]-1)*(GammaEff[i]+1)*dmdR*c**2) < (GammaEff[i]*(gammaAdi[i]-1)*Eint2[i]/V2[i]*(dmdR/rhoprim[i] + m[i]/rhoprim[i]*dlnrhodR)):
                    print '\n\ndmdR=%s\nGamma[i]=%s\nGammaEff[i]=%s\ns/R[i]=%s\ndmdR\nm[i]=%s\ni=%d\nGamma0=%s\nrho=%s'%(dmdR,Gamma[i],GammaEff[i],s/R[i],dmdR/m[i],i
,Gamma0,rho[i])

                                         

#                #Eint[i] = -4/3.*R[i]**3*pi*c**2*dGammaEffdGamma[i] * (M0+m[i]) * (mup * ((Gamma[i]-1)*epsilonp/mup + 1) + mue * ((Gamma[i]-1)*epsilone/mue + 1))
#                else:
#                    dEaddRterm = -GammaEff[i] * np.sum(Rsum2*rho[:i+1]*Rfraq * Gammafraq*(mup*pp2 / RGFraqSqrtP + mue*pe2/RGFraqSqrtE) * (np.cos(theta[:i+1])-1) * Rdist) #Note from 16/8
#                    Bterm = GammaEff[i] * np.sum(rho[:i+1]*(mup*RGFraqSqrtP + mue*RGFraqSqrtE)*np.sin(theta[:i+1])*dthetadR*Rsum2*Rdist)

#                    Gamma[i+1] = Gamma[i] - ((GammaEff[i]+1)*(Gamma[i]-1)*c**2*dmdR - (2*pi*c**2/R[i] + np.sin(theta[i])/(1-np.cos(theta[i]))*dthetadR[i] )*dEaddRterm - 2*pi*c**2*Bterm) / ((M0+m[i])*c**2 + Eint2[i] + 4*pi*c**2/3./Gamma[i] * dEaddRterm) * (R[i+1] - R[i]) #Note from 25/8-13

#                    dEaddR[i] = np.copy(dEaddRterm)
#                    Bterm_array[i] = np.copy(Bterm)

        if A0==0:print Gamma[i+1]
       
        gammaAdi[i+1] = (4+1/Gamma[i+1]) / 3.

        GammaEff[i+1] = (gammaAdi[i+1]*Gamma[i+1]**2 - gammaAdi[i+1]+1)/Gamma[i+1]

        dGammaEffdGamma[i+1] = (4+1/Gamma[i+1]**2 + 2/Gamma[i+1]**3) / 3.

        beta[i+1] = (1-1/Gamma[i+1]**2)**(1/2.)
        

        if reverseShock:
            if Gamma[0]/Gamma[i+1] > 1.001:  
                bbterm = beta0 * beta[i+1]
                gamma43[i+1] = Gamma0*Gamma[i+1]*(1 - bbterm)
                if gamma43[i+1] <= 1.: gamma43[i+1] = 1+1e-5
            else: gamma43[i+1] = 1.+1e-5
            gammaAdi3[i+1] = (4+1/gamma43[i+1]) / 3.    #Is it supposed to be the same as for the FS?
            beta43[i+1] = (1-1/gamma43[i+1]**2)**0.5
            GammaEff3[i+1] = (gammaAdi3[i+1]*Gamma[i+1]**2 - gammaAdi3[i+1]+1)/Gamma[i+1]
            dGammaEff3dGamma[i+1] = (4+1/Gamma[i+1]**2 + 2/Gamma[i+1]**3) / 3.
        if Gamma[i+1]<1: 
            return grandNan

        tburst[i+1] = tburst[i] + (R[i+1]-R[i]) / beta[i]/c     #Burster time adder
        tcomoving[i+1] = tcomoving[i] + (tburst[i+1]-tburst[i]) / Gamma[i]         #Comoving time adder
        if A0 == 0: tobs[i+1] = tobs[i] + (tburst[i+1] - tburst[i]) * (1+z) / ((1+beta[i])*Gamma[i]**2)
        else:tobs[i+1] = tobs[i] + (tburst[i+1] - tburst[i]) * (1+z) / ((1+beta[i])*Gamma[i]**2) ### (1+z)*(tburst[i+1] - (R[i+1])/c)   ###  ###
        if useSpread:# and (R[i] > startSpread): 
            if theta[i] >= np.pi / 2: theta[i+1] = np.pi/2
            else:
#                gammabar = (4*Gamma[i]+1) / (3*Gamma[i])
#                vperp = ( gammabar*(gammabar-1)*(Gamma[i]-1) / (1 + gammabar*(Gamma[i]-1)) )**(1/2.) * c
                theta[i+1] = theta[i] + dthetadR[i] * (R[i+1] - R[i])    ###2*vperp / (R[i]+R[i+1]) * (tcomoving[i+1]-tcomoving[i])    #Spreading jet



        else: theta[i+1] = theta[i]
        if A0==0: #Cocoon
            r0 = 1e7
                                         

#            for i_coc in range(1,i+1):   ### Calculating the density of the cocoon for each radial step. Only used when photosphere==True                  
            
            rho_coc[i] = M0 / (2 * np.pi * np.sum((R[1:i+(i==0)+1]-R[0])**2*Rcoc_dist[1:i+(i==0)+1]*(1-np.cos(theta[1:i+(i==0)+1])))) ### Factor i+(i==0) is to set index 0 to the same value as index 1. Otherwise we have a singularity at 0
                #


                #Setting these two first lines to see evaluation                                                                                                        
            kappa_0 =  5e24
            kappa_R = kappa_0 * rho_coc[i] * T[i]**(-7/2.)
#            print R
#            print rho_coc[i]
            Dtau = tau_coc / rho_coc[i] / kappa_R   ### The distance with optical depth tau=1                                                                          
            print R[i]
            print 'Gamma = %s'%Gamma[i]
            print 'Gamma[i+1] = %s'%Gamma[i+1]
            R_photo[i] = R[i] - Dtau    ### The radius of the photosphere
#            print Dtau
#            print R_photo[i]
            if True: ### Description from Pe'er at al 2007
                T[i] = (E0 / 4 / np.pi / r0**2 / c / rad_const) ** 0.25 / Gamma[i] * (R_photo[i] / r0) ** (-2/3.)
                print T[i]

            Acoc = 2*np.pi * R[i] ** 2 * (1-np.cos(theta[i]))   #Area of cocoon
            print gammap[i]
            print pe2[i]
            print gammae[i]
            print Gamma[i]
            print dGammadR
            dEaddR = -c**2 * (1/R[i] + 1/3/Gamma[i] * dGammadR) * M0 * ( mup * pp2[i] / gammap[i] + mue * pe2[i] / gammae[i])
            print -sigmaT * Acoc * T[i]**4 / Gamma[i] / beta[i] / c
            raw_input(dEaddR)
            Eint2[i+1] = Eint2[i] + ( -sigmaT * Acoc * T[i]**4 / Gamma[i] / beta[i] / c  +  dEaddR ) * (R[i+1] - R[i])
            
##            T[i+1] = T[i] - sigma_B * T[i]**4 * 2*np.pi*(1-np.cos(2*theta[i]))*R[i]**2 / M0 * (me+mp) / 3 / kB * (tcomoving[i+1]-tcomoving[i])
#                 else:
 #                   T[i+1] = Eint[i+1] / M0 * (me+mp) / 3 / kB
#            print 'T = %s'%T[i+1]
#            print 'tobs = %s'%tobs[i+1]
#            raw_input('R = %s'%R[i+1])
#            T[i+1] = T[i] + dEaddR / M0 * (me+mp) / 3 / kB * (R[i+1] - R[i])
#            T[i+1] = np.copy(T[i])


            gammae[i+1] = Eint2[i+1] / 2 / M0 / mue / c**2 + 1
            gammap[i+1] = mue / mup * (gammae[i+1] - 1) + 1

            if (gammae[i+1] < 1) or (gammap[i+1] < 1): 
                return grandNan

            pp2[i+1] = np.sqrt(gammap[i+1]**2 - 1)
            pe2[i+1] = np.sqrt(gammae[i+1]**2 - 1)

            
        #Injection of mass. Note from 11/2 2014
        if reverseShock:
            if ( M3[i] <= M0 ):  
                dM3 = M0 * (beta0 - beta[i]) / beta[i] / tprompt * (tcomoving[i+1] - tcomoving[i])
                rho4[i+1] = M0 / (2*pi*R[i+1]**2*Delta0*(1-np.cos(theta[i+1]) )*beta[0]*c)
            else: 
                if shutOff == False: 
                    RScutOff = tcomoving[i]
                    if printStats: print "Reverse shock feed stopped after %s s, observer\'s time"%tobs[i]
                    shutOff = True
                    dM3 = M0 - M3[i]
                dM3 = 0.
        
            M3[i+1] = M3[i] + dM3
            dM3dR = dM3/(R[i+1] - R[i])
        
        #Check if we have passed the latest data entry time
        if tobs[i-3] > tobsEnd:     #Even if i == 0, tobs[i-1]=tobs[-1] = 0 from space allocation
            #Cleaning...
            Gamma = Gamma[:i+1]
#            Eint2 = Eint2[:i+1]
            tburst = tburst[:i+1]
            tobs = tobs[:i+1]
            tcomoving = tcomoving[:i+1]
            theta = theta[:i+1]
            m = m[:i+1]
            if A0==0:
                T = T[:i+1]
            R = R[:i+1]
            gammae = gammae[:i+1]
            if reverseShock:
                gamma43 = gamma43[:i+1]
                M3 = M3[:i+1]
                rho4 = rho4[:i+1]
                gammaeRS = gammaeRS[:i+1]
            
            break
        
        if i == (len(R)-2): #We need more elements
            R = np.append(R,R[-1]**2/R[-2])
            Gamma = np.append(Gamma,0.)
            Eint2 = np.append(Eint2,0.)
            tburst = np.append(tburst,0.)
            tobs = np.append(tobs,0.)
            tcomoving = np.append(tcomoving,0.)
            theta = np.append(theta,0.)
            m = np.append(m,0.)
            beta = np.append(beta,0.)
#            dEaddR = np.append(dEaddR,0.)
#            Bterm_array = np.append(dEaddR,0.)
            
            gammaAdi = np.append(gammaAdi,0.)
            GammaEff = np.append(GammaEff,0.)
            dGammaEffdGamma = np.append(dGammaEffdGamma,0.)
            gammae = np.append(gammae,0.)
            gammap = np.append(gammap,0.)
            pp2 = np.append(pp2,0.)
            pe2 = np.append(pe2,0.)
            rhoprim = np.append(rhoprim,0.)
            V2 = np.append(V2,0.)
            epsilon = np.append(epsilon,0.)
            if reverseShock:
                V3 = np.append(V3,0.)
                gamma43 = np.append(gamma43,0.)
                beta43 = np.append(beta43,0.)
                dGammaEff3dGamma = np.append(dGammaEff3dGamma,0.)
                gammaAdi3 = np.append(gammaAdi3,0.)
                GammaEff3 = np.append(GammaEff3,0.)
                M3 = np.append(M3,0.)
                rho4 = np.append(rho4,0.)
                Eint3 = np.append(Eint3,0.)
                gammaeRS = np.append(gammaeRS,0.)
                gammapRS = np.append(gammapRS,0.)
                pp2rs = np.append(pp2rs,0.)
                pe2rs = np.append(pe2rs,0.)
                epsilonRS = np.append(epsilonRS,0.)


            if (s != 0) and (R[i] > 10**profile_cutoff): rho = np.append(rho,A0*10**(-profile_cutoff*s))
            else:rho = np.append(rho,A0 * R[-1]**(-s))
            if A0==0:
                T = np.append(T,0.)
        i += 1
    if RScutOff == 0: RScutOff = tcomoving[-1]
    if reverseShock == False: RScutOff = 0. 
    if A0==0: #Cocoon
        print i
        print Gamma
        return Gamma,tburst,tobs,tcomoving,theta,M0+m,m,M3,rho4,T,R
    else: #Jet
#        plt.plot(tobs,Eint)
#        plt.loglog()
#        plt.show()

#        print np.shape(R)
#        print np.shape(dthetadR)
#        plt.plot(tobs[1:],dthetadR[1:]*(R[1:]-R[:-1]))
#        plt.loglog()
#        plt.show()

#        plt.subplot(1,2,1)
#        plt.plot(tobs,dEaddR[:len(tburst)])
#        plt.loglog()

#        plt.subplot(1,2,2)

#        plt.plot(tobs,Bterm_array[:len(tburst)])
#        plt.loglog()
#        plt.show()

#        plt.plot(tobs,theta)
#        plt.xscale('log')
#        plt.show()

#        plt.plot(R,Gamma)
#        plt.loglog()
#        plt.show()

#        np.savetxt('gamma_old.txt',[tobs,Gamma])

        if reverseShock: return Gamma,tburst,tobs,tcomoving,theta,M0+m,m,M3,rho4,RScutOff,R,gamma43,gammae,gammaeRS
        else: return Gamma,tburst,tobs,tcomoving,theta,M0+m,m,R,gammae
 
