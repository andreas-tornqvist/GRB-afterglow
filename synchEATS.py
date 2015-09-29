def specGenEATS(Gamma,GammaEff,t,eB,eE,E52,nprim,p,nu,B,kappas,opticalDepth,depth):
    #The GammaEff is the internal electron lorentz factor. This should depend on the bulk lorentz factor
    import numpy as np
    c,pi,mp,me,qe, sigmaT = 2.9979e10,np.pi, 1.6726e-24, 9.1094e-28, 4.803204e-10, 6.6524e-25
#This module considers synchrotron emission in the fluid frame according to the Ph.D. Thesis of G. Johannesson
    
    #phipF,phipS,XpF,XpS,kappa1,kappa2,kappa3,kappa4 = opDepthConst    #Constants to optical depth calculation

    #nprim = n * 4 * Gamma    #Shock front density 

    #Pmax = me*c**2*sigmaT*Gamma*B/3/qe
    #phipF = 1.89 - 0.935*p + 0.17*p**2
    #phipS = 0.54 + 0.08*p
    #kappas = [-kappa1/3,kappa1/2,-1/kappa1,kappa2*(p-1)/2,-1/kappa2,-kappa3/3,kappa3*(p-1)/2,-1/kappa3,kappa4/2,-1/kappa4]
    #if True:
    #PmaxF = phipF*2.234*qe**3*nprim*B/me/c**2  
    PmaxF = kappas[10]*2.234*qe**3*nprim*B/me/(c**2)
    #PmaxF = kappas[10]*2.234*qe**3*B/me/c**2
    #print n
    #raw_input(PmaxF)
    PmaxS = kappas[11]*11.17*(p-1)*qe**3*nprim*B/(3*p-1)/me/(c**2)
    #PmaxS = kappas[11]*11.17*(p-1)*qe**3*B/(3*p-1)/me/c**2
    #raw_input(PmaxS)


    #else:
    #    PmaxF = me * c**2 * sigmaT * B / 3 / qe / Gamma  *  nprim
    #    PmaxS = PmaxF

    XpF = 0.455 + 0.08*p
    XpS = 0.06 + 0.28*p


    #Minimum and characteristic frequencies

    #gammam = (p-2)/(p-1) * (eE*mp/me*(Gamma-1)+1)
    #gammac = 6*pi*me*c/sigmaT/Gamma/B**2/t

    #num = Gamma * gammam**2 *qe * B / 2 / pi / me / c
    #nuc = Gamma * gammac**2 *qe * B / 2 / pi / me / c
    #print t
    
    gammamFac1 = (p-2)/(p-1)
    gammamFac2 = eE*mp/me
    gammam = gammamFac1 * (gammamFac2 * (GammaEff -1) + 1)
    #print gammam
    #print Gamma
    #print B
    #print t
    #raw_input("...")
    #gammac = 6*pi*me*c/sigmaT/Gamma/B**2/t
    pimec = 2 * pi * me * c
    gammac1 = 3*pimec/sigmaT
    B2t = B**2*t
    gammac = gammac1/Gamma/B2t
    #raw_input(gammac)
    #print GammaEff
    #print gammam
    #print gammac

    gamToNuFactor = qe * B / pimec
    num = gammam**2 * gamToNuFactor
    nuc = gammac**2 * gamToNuFactor
    
    #num = (eE*mp/me * (GammaEff-1) + 1) ** 2 * ((p-2)/(p-1))**2            #Typical frequency, comoving frame
    #nuc = 18.*pi*qe*me*c/(sigmaT**2*B**3*t**2)                             #Cooling frequency, comoving frame
    
    #print len(nu)
    #print len(Gamma)
    #print len(num)
    #print len(nuc)

    #The Sari, Piran & Narayan (1998) approach:
    #P = PmaxF * (num>nuc)*((nuc>nu)*(nu/nuc)**(1/3.) + ((num>nu)&(nu>=nuc))*(nu/nuc)**(-1/2.) + (nu>=num)*(num/nuc)**(-1/2.)*(nu/num)**(-p/2.) )    +     PmaxS(num<=nuc)*((num>nu)*(nu/num)**(1/3.) + ((nuc>nu)&(nu>=num))*(nu/num)**(-(p-1)/2.) + (nu>nuc)*(nu/num)**(-(p-1)/2.)*(nu/nuc)**(-p/2.))
    #kappa1,kappa2,kappa3,kappa4 = 2.37-0.3*p  ,  14.7 - 8.68*p + 1.4*p**2  ,  6.94 - 3.844*p + 0.62*p**2  ,  3.5 - 0.2*p
    if False: #calculate spectrum in for-loop to gain speed. Does not work; test shows that the array multiplyer works faster.
        P = np.zeros(len(nu))
        numnuc = num > nuc
        nunum = nu/num
        nunuc = nu/nuc
        #kappas[0] = -kappa1/3
        #kappa[1] = kappa1/2
        #kappas[2] = -1/kappa1
        #kappa2p = kappa2*(p-1)/2
        #kappas[4] = -1/kappa2
        #kappas[5] = -kappa3/3
        #kappa[6] = kappa3*(p-1)/2
        #kappas[7] = -1/kappa3
        #kappas[8] = kappa4/2
        #kappas[9] = -1/kappa4
        for nui in range(len(nu)):
            if numnuc[nui]:
                P[nui] = PmaxF[nui] * (nunuc[nui]**kappas[0] + nunuc[nui]**kappas[1]) ** kappas[2] * (1+nunum[nui]**kappas[3])**kappas[4]
            else:
                P[nui] = PmaxS[nui] * (nunum[nui]**kappas[5] + nunum[nui]**kappas[6])**kappas[7] * (1+nunuc[nui]**kappas[8])**kappas[9]

    #The Jr approach
    
    #print np.size((num>nuc) * ((nu/nuc)**(-kappa1/3)))
    #print np.size((nu/nuc)**(kappa1/2) ** (-1/kappa1) * (1+(nu/num)**(kappa2*(p-1)/2))**(-1/kappa2))
    #print np.size((num<=nuc) * ((nu/num)**(-kappa3/3)))
    #print np.size((nu/num)**(kappa3*(p-1)/2))**(-1/kappa3) * (1+(nu/nuc)**(kappa4/2))**(-1/kappa4))
    else:
        #print kappa1
        #print kappa2
        #print kappa3
        #print kappa4
        #print '...---...'
        #print nu
        #print num
        #print nuc
        #print (num>nuc) * ((nu/nuc)**(-kappa1/3) + (nu/nuc)**(kappa1/2)) ** (-1/kappa1) 
        #print (1+(nu/num)**(kappa2*(p-1)/2))**(-1/kappa2) 
        #raw_input((num<=nuc) * ((nu/num)**(-kappa3/3) + (nu/num)**(kappa3*(p-1)/2))**(-1/kappa3) * (1+(nu/nuc)**(kappa4/2))**(-1/kappa4))
        #term11 = num>nuc
        #nunuc = nu/nuc
        #nunum = nu/num
        #term12 = nunuc**(kappas[0])
        #term21 = nunuc**(kappas[1])
        #term31 = 1+nunum**(kappas[3])
        #term41 = nunum**(kappas[5])
        #term42 = nunum**(kappas[6])
        #term51 = (1+(nu/nuc)**(kappas[8]))
        #P = np.zeros(len(nu))
        #PmaxFRed = PmaxF * term11
        #PmaxSRed = PmaxS * (term11==False)
        #P[term11] = PmaxF * (term12 + term21) ** (kappas[2]) * term31 ** (kappas[4])
        #P[term11==False] = PmaxS * (term41 + term42) ** (kappas[7]) * term51 ** (kappas[9])
        P = PmaxF * (num>nuc) * ((nu/nuc)**(kappas[0]) + (nu/nuc)**(kappas[1])) ** (kappas[2]) * (1+(nu/num)**(kappas[3]))**(kappas[4])   +   PmaxS * (num<=nuc) * ((nu/num)**(kappas[5]) + (nu/num)**(kappas[6]))**(kappas[7]) * (1+(nu/nuc)**(kappas[8]))**(kappas[9])
        #P = PmaxF * (num>nuc) * ((nu/nuc)**(-kappa1/3) + (nu/nuc)**(kappa1/2)) ** (-1/kappa1) * (1+(nu/num)**(kappa2*(p-1)/2))**(-1/kappa2)   +   PmaxS * (num<=nuc) * ((nu/num)**(-kappa3/3) + (nu/num)**(kappa3*(p-1)/2))**(-1/kappa3) * (1+(nu/nuc)**(kappa4/2))**(-1/kappa4)
#F = Fmax * ((num>nuc)*((nuc>nu)*(nu/nuc)**(1/3.) + ((num>nu)&(nu>=nuc))*(nu/nuc)**(-1/2.) + (nu>=num)*(num/nuc)**(-1/2.)*(nu/num)**(-p/2.) ) + (num<=nuc)*((num>nu)*(nu/num)**(1/3.) + ((nuc>nu)&(nu>=num))*(nu/num)**(-(p-1)/2.) + (nu>nuc)*(nu/num)**(-(p-1)/2.)*(nu/nuc)**(-p/2.)))
    #print np.shape(P)


    

    ######################
    #   Optical depth    #
    ######################

    if opticalDepth:
        #PmaxF = phipF*2.234*qe**3*nprim*B/me/c**2  
        #PmaxS = phipS*11.17*(p-1)*qe**3*nprim*B/(3*p-1)/me/c**2
    
        alpha0F = 11.7 * kappas[10] * XpF**(-3) * qe * nprim /  B  *  gammac**(-5)
        alpha0S =  7.8 * kappas[11] * XpS**(-(4+p)/2.) * (p+2)*(p-1)* qe * nprim * gammam**(-5) / (p+2/3.)/  B
        alphanu = alpha0F * (num>nuc) * ( (nu<=nuc) * (nu/nuc)**(-5/3.) + ((nuc<nu)&(nu<num)) * (nu/nuc)**(-3) + (num<=nu) * (num/nuc)**-3*(nu/num)**(-(p+5)/2.) )    +    alpha0S * (num<=nuc) * ( (nu<=num)*(nu/num)**(-5/3.) + ((num<nu)&(nu<nuc)) * (nu/num)**(-(p+4)/2.) + (nuc<=nu) * (nuc/num)**(-(p+4)/2.) * (nu/nuc)**(-(p+5)/2.) )

        tau = alphanu * depth     #tau is always negative
#        print tau
        tau[tau>-1e-8] = 0.
        P = P * np.exp(tau)


    
    return P
