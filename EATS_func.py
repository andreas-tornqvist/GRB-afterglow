def eats_function(p , nu, slow_cooling,fast_cooling, num, nuc, PmaxF, PmaxS):
    import numpy as np

 #Constants to optical depth calculation
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

    P_out = np.zeros(len(PmaxF))

#    P_out[fast_cooling] = PmaxF[fast_cooling] * ((nu[fast_cooling]/nuc[fast_cooling])**(kappa13) + (nu[fast_cooling]/nuc[fast_cooling])**(kappa12)) ** (kappa11) * (1+(nu[fast_cooling]/num[fast_cooling])**(kappa2p))**(kappa12inv)

#    P_out[slow_cooling] = PmaxS[slow_cooling] * ((nu[slow_cooling]/num[slow_cooling])**(kappa33) + (nu[slow_cooling]/num[slow_cooling])**(kappa3p))**(kappa13inv) * (1+(nu[slow_cooling]/nuc[slow_cooling])**(kappa42))**(kappa14)


    P_out[fast_cooling] = PmaxF[fast_cooling] * ((nu[fast_cooling]/nuc[fast_cooling])**(kappa13) + (nu[fast_cooling]/nuc[fast_cooling])**(kappa12)) ** (kappa11) * (1+(nu[fast_cooling]/num[fast_cooling])**(kappa2p))**(kappa12inv)

    P_out[slow_cooling] =  PmaxS[slow_cooling] * ((nu[slow_cooling]/num[slow_cooling])**(kappa33) + (nu[slow_cooling]/num[slow_cooling])**(kappa3p))**(kappa13inv) * (1+(nu[slow_cooling]/nuc[slow_cooling])**(kappa42))**(kappa14)


#    print P_out
                
    return P_out
