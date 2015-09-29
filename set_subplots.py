def set_subplots(reverseShock,n_params,parametrar):
    import numpy as np
    from matplotlib import gridspec
    from matplotlib import pyplot as plt

    #Plot labels                                                                                                                                                        
    epsilone_label = r'$\epsilon_{\rm e}$'
    epsilone_FS_label = r'$\epsilon_{\rm e,FS}$'
    epsilone_RS_label = r'$\epsilon_{\rm e,RS}$'
    E0_label = r'$E_0[{\rm erg}]$'
    n_label = r'$n[{\rm cm}^{-3}]$'
    Gamma0_label = r'$\Gamma_0$'
    epsilonB_label = r'$\epsilon_{\rm B}$'
    epsilonB_FS_label = r'$\epsilon_{\rm B,FS}$'
    epsilonB_RS_label = r'$\epsilon_{\rm B,RS}$'
    p_label = r'$p$'
    p_FS_label = r'$p_{\rm FS}$'
    p_RS_label = r'$p_{\rm RS}$'
    theta0_label = r'$\theta_0$'
    alpha_label = r'$\alpha$'
    tprompt_label = r'$\Delta t_{\rm prompt}[s]$'
    T_coc_label = r'$T_{\rm coc}$ [K]'
    R_coc_label = r'$R_{\rm coc}$ [cm]'
    theta0_coc_label = r'$\theta_{\rm 0,coc}$ [$^{\circ}$]'
    Gamma0_coc_label = r'$\Gamma_{\rm 0,coc}$'
    N_coc_label = r'$N_{\rm coc}$'

    if not reverseShock:
        if n_params <= 8:
            plt.figure(figsize=(5*n_params, 2*n_params))
            G = gridspec.GridSpec(2,4)
            epsilone_subplot = plt.subplot(G[0,3],xscale=u'log',xlabel=epsilone_label,ylabel='Probability')
            E0_subplot = plt.subplot(G[0,0],xscale=u'log',xlabel=E0_label,ylabel='Probability')
            n_subplot = plt.subplot(G[0,1],xscale=u'log',xlabel=n_label,ylabel='Probability')
            Gamma0_subplot = plt.subplot(G[1,1],xscale=u'log',xlabel=Gamma0_label,ylabel='Probability')
            epsilonB_subplot = plt.subplot(G[1,3],xscale=u'log',xlabel=epsilonB_label,ylabel='Probability')
            p_subplot = plt.subplot(G[1,0],xlabel=p_label,ylabel='Probability')
            theta0_subplot = plt.subplot(G[0,2],xlabel=theta0_label,ylabel='Probability')
            alpha_subplot = plt.subplot(G[1,2],xlabel=alpha_label,ylabel='Probability')
            T_coc_subplot = None
            R_coc_subplot = None
            theta0_coc_subplot = None
            Gamma0_coc_subplot = None
            N_coc_subplot = None

        elif n_params >= 11: #4 Thermal parameters
            plt.figure(figsize=(5*n_params, 2*n_params))
            G = gridspec.GridSpec(3,4)
            epsilone_subplot = plt.subplot(G[0,3],xscale=u'log',xlabel=epsilone_label,ylabel='Probability')
            E0_subplot = plt.subplot(G[0,0],xscale=u'log',xlabel=E0_label,ylabel='Probability')
            n_subplot = plt.subplot(G[0,1],xscale=u'log',xlabel=n_label,ylabel='Probability')
            Gamma0_subplot = plt.subplot(G[1,1],xscale=u'log',xlabel=Gamma0_label,ylabel='Probability')
            epsilonB_subplot = plt.subplot(G[1,3],xscale=u'log',xlabel=epsilonB_label,ylabel='Probability')
            p_subplot = plt.subplot(G[1,0],xlabel=p_label,ylabel='Probability')
            theta0_subplot = plt.subplot(G[0,2],xlabel=theta0_label,ylabel='Probability')
            alpha_subplot = plt.subplot(G[1,2],xlabel=alpha_label,ylabel='Probability')
            T_coc_subplot = plt.subplot(G[2,0],xscale=u'log',xlabel=T_coc_label,ylabel='Probability')
            R_coc_subplot = plt.subplot(G[2,1],xscale=u'log',xlabel=R_coc_label,ylabel='Probability')
            theta0_coc_subplot = plt.subplot(G[2,2],xscale=u'log',xlabel=theta0_coc_label,ylabel='Probability')
            Gamma0_coc_subplot = plt.subplot(G[2,3],xscale=u'log',xlabel=Gamma0_coc_label,ylabel='Probability')
            N_coc_subplot = plt.subplot(G[2,2],xscale=u'log',xlabel=N_coc_label,ylabel='Probability')
        

        epsilone_subplot.xaxis.label.set_fontsize(50)
        epsilone_subplot.yaxis.label.set_fontsize(40)
        E0_subplot.xaxis.label.set_fontsize(50)
        E0_subplot.yaxis.label.set_fontsize(40)
        n_subplot.xaxis.label.set_fontsize(50)
        n_subplot.yaxis.label.set_fontsize(40)
        Gamma0_subplot.xaxis.label.set_fontsize(50)
        Gamma0_subplot.yaxis.label.set_fontsize(40)
        epsilonB_subplot.xaxis.label.set_fontsize(50)
        epsilonB_subplot.yaxis.label.set_fontsize(40)
        p_subplot.xaxis.label.set_fontsize(50)
        p_subplot.yaxis.label.set_fontsize(40)
        theta0_subplot.xaxis.label.set_fontsize(50)
        theta0_subplot.yaxis.label.set_fontsize(40)
        s_subplot = None
        if n_params >= 11:
            T_coc_subplot.xaxis.label.set_fontsize(50)
            T_coc_subplot.yaxis.label.set_fontsize(50)
            R_coc_subplot.xaxis.label.set_fontsize(50)
            R_coc_subplot.yaxis.label.set_fontsize(50)
            theta0_coc_subplot.xaxis.label.set_fontsize(50)
            theta0_coc_subplot.yaxis.label.set_fontsize(50)
            Gamma0_coc_subplot.xaxis.label.set_fontsize(50)
            Gamma0_coc_subplot.yaxis.label.set_fontsize(50)


        alpha_subplot.xaxis.label.set_fontsize(50)
        alpha_subplot.yaxis.label.set_fontsize(40)

        epsilone_RS_subplot = None
        epsilonB_RS_subplot = None
        tprompt_subplot = None
        p_RS_subplot = None
#        return epsilone_subplot, E0_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, alpha_subplot
    else: 
        if n_params == 11: #Horisontal orientaion, 11 parameters
            plt.figure(figsize=(5*n_params, 2*n_params))
            G = gridspec.GridSpec(6,5)
            epsilone_subplot = plt.subplot(G[:3,3],xscale=u'log',xlabel=epsilone_FS_label,ylabel='Probability')
            epsilone_RS_subplot = plt.subplot(G[3:,3],xscale=u'log',xlabel=epsilone_RS_label,ylabel='Probability')
            E0_subplot = plt.subplot(G[:3,0],xscale=u'log',xlabel=E0_label,ylabel='Probability')
            n_subplot = plt.subplot(G[:3,1],xscale=u'log',xlabel=n_label,ylabel='Probability')
            Gamma0_subplot = plt.subplot(G[3:,0],xscale=u'log',xlabel=Gamma0_label,ylabel='Probability')
            epsilonB_subplot = plt.subplot(G[:2,4],xscale=u'log',xlabel=epsilonB_FS_label,ylabel='Probability')
            epsilonB_RS_subplot = plt.subplot(G[2:4,4],xscale=u'log',xlabel=epsilonB_RS_label,ylabel='Probability')
            p_subplot = plt.subplot(G[3:,1],xlabel=p_FS_label,ylabel='Probability')
            p_RS_subplot = plt.subplot(G[3:,2],xlabel=p_RS_label,ylabel='Probability')
            theta0_subplot = plt.subplot(G[:3,2],xlabel=theta0_label,ylabel='Probability')
            tprompt_subplot = plt.subplot(G[4:,4],xlabel=tprompt_label,ylabel='Probability')

            epsilone_subplot.xaxis.label.set_fontsize(50)
            epsilone_subplot.yaxis.label.set_fontsize(40)
            epsilone_RS_subplot.xaxis.label.set_fontsize(50)
            epsilone_RS_subplot.yaxis.label.set_fontsize(40)
            E0_subplot.xaxis.label.set_fontsize(50)
            E0_subplot.yaxis.label.set_fontsize(40)
            n_subplot.xaxis.label.set_fontsize(50)
            n_subplot.yaxis.label.set_fontsize(40)
            Gamma0_subplot.xaxis.label.set_fontsize(50)
            Gamma0_subplot.yaxis.label.set_fontsize(40)
            epsilonB_subplot.xaxis.label.set_fontsize(50)
            epsilonB_subplot.yaxis.label.set_fontsize(40)
            p_subplot.xaxis.label.set_fontsize(50)
            p_subplot.yaxis.label.set_fontsize(40)
            epsilonB_RS_subplot.xaxis.label.set_fontsize(50)
            epsilonB_RS_subplot.yaxis.label.set_fontsize(40)
            p_RS_subplot.xaxis.label.set_fontsize(50)
            p_RS_subplot.yaxis.label.set_fontsize(40)
            theta0_subplot.xaxis.label.set_fontsize(50)
            theta0_subplot.yaxis.label.set_fontsize(40)
            tprompt_subplot.xaxis.label.set_fontsize(50)
            tprompt_subplot.yaxis.label.set_fontsize(40)
            alpha_subplot = None
            s_subplot = None
            N_coc_subplot = None
            T_coc_subplot = None
            R_coc_subplot = None
            theta0_coc_subplot = None
            Gamma0_coc_subplot = None

#            return epsilone_subplot, E0_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, epsilonB_RS_subplot, p_RS_subplot, epsilone_RS_subplot, tprompt_subplot

        elif n_params == 10: #Horisontal orientation, 10 parameters
            plt.figure(figsize=(5*n_params, 2*n_params))
            G = gridspec.GridSpec(2,5)
            epsilone_subplot = plt.subplot(G[0,3],xscale=u'log',xlabel=epsilone_FS_label,ylabel='Probability')
            epsilone_RS_subplot = plt.subplot(G[1,3],xscale=u'log',xlabel=epsilone_RS_label,ylabel='Probability')
#            E0_subplot = plt.subplot(G[:3,0],xscale=u'log',xlabel=E0_label,ylabel='Probability')
            n_subplot = plt.subplot(G[0,0],xscale=u'log',xlabel=n_label,ylabel='Probability')
            Gamma0_subplot = plt.subplot(G[0,1],xscale=u'log',xlabel=Gamma0_label,ylabel='Probability')
            epsilonB_subplot = plt.subplot(G[0,4],xscale=u'log',xlabel=epsilonB_FS_label,ylabel='Probability')
            epsilonB_RS_subplot = plt.subplot(G[1,4],xscale=u'log',xlabel=epsilonB_RS_label,ylabel='Probability')
            p_subplot = plt.subplot(G[1,0],xlabel=p_FS_label,ylabel='Probability')
            p_RS_subplot = plt.subplot(G[1,1],xlabel=p_RS_label,ylabel='Probability')
            theta0_subplot = plt.subplot(G[1,2],xlabel=theta0_label,ylabel='Probability')
            tprompt_subplot = plt.subplot(G[0,2],xlabel=tprompt_label,xscale=u'log',ylabel='Probability')
            
            epsilone_subplot.xaxis.label.set_fontsize(50)
            epsilone_subplot.yaxis.label.set_fontsize(40)
            epsilone_RS_subplot.xaxis.label.set_fontsize(50)
            epsilone_RS_subplot.yaxis.label.set_fontsize(40)
            n_subplot.xaxis.label.set_fontsize(50)
            n_subplot.yaxis.label.set_fontsize(40)
            Gamma0_subplot.xaxis.label.set_fontsize(50)
            Gamma0_subplot.yaxis.label.set_fontsize(40)
            epsilonB_subplot.xaxis.label.set_fontsize(50)
            epsilonB_subplot.yaxis.label.set_fontsize(40)
            p_subplot.xaxis.label.set_fontsize(50)
            p_subplot.yaxis.label.set_fontsize(40)
            epsilonB_RS_subplot.xaxis.label.set_fontsize(50)
            epsilonB_RS_subplot.yaxis.label.set_fontsize(40)
            p_RS_subplot.xaxis.label.set_fontsize(50)
            p_RS_subplot.yaxis.label.set_fontsize(40)
            theta0_subplot.xaxis.label.set_fontsize(50)
            theta0_subplot.yaxis.label.set_fontsize(40)
            tprompt_subplot.xaxis.label.set_fontsize(50)
            tprompt_subplot.yaxis.label.set_fontsize(40)
            s_subplot = None
            N_coc_subplot = None
            T_coc_subplot = None
            R_coc_subplot = None
            theta0_coc_subplot = None
            Gamma0_coc_subplot = None
            E0_subplot = None
            alpha_subplot = None

#            return epsilone_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, alpha_subplot, epsilonB_RS_subplot, p_RS_subplot, epsilone_RS_subplot, tprompt_subplot
        elif n_params==12:
            plt.figure(figsize=(5*n_params, 2*n_params))
            G = gridspec.GridSpec(3,4)
            epsilone_subplot = plt.subplot(G[0,3],xscale=u'log',xlabel=epsilone_FS_label,ylabel='Probability')
            epsilone_RS_subplot = plt.subplot(G[1,3],xscale=u'log',xlabel=epsilone_RS_label,ylabel='Probability')
            E0_subplot = plt.subplot(G[0,0],xscale=u'log',xlabel=E0_label,ylabel='Probability')
            n_subplot = plt.subplot(G[0,1],xscale=u'log',xlabel=n_label,ylabel='Probability')
            Gamma0_subplot = plt.subplot(G[1,0],xscale=u'log',xlabel=Gamma0_label,ylabel='Probability')
            epsilonB_subplot = plt.subplot(G[0,2],xscale=u'log',xlabel=epsilonB_FS_label,ylabel='Probability')
            epsilonB_RS_subplot = plt.subplot(G[1,2],xscale=u'log',xlabel=epsilonB_RS_label,ylabel='Probability')
            p_subplot = plt.subplot(G[2,0],xlabel=p_FS_label,ylabel='Probability')
            p_RS_subplot = plt.subplot(G[2,1],xlabel=p_RS_label,ylabel='Probability')
            theta0_subplot = plt.subplot(G[2,2],xlabel=theta0_label,ylabel='Probability')
            alpha_subplot = plt.subplot(G[2,3],xlabel=alpha_label,ylabel='Probability')
            tprompt_subplot = plt.subplot(G[1,1],xlabel=tprompt_label,ylabel='Probability')

            epsilone_subplot.xaxis.label.set_fontsize(50)
            epsilone_subplot.yaxis.label.set_fontsize(40)
            epsilone_RS_subplot.xaxis.label.set_fontsize(50)
            epsilone_RS_subplot.yaxis.label.set_fontsize(40)
            E0_subplot.xaxis.label.set_fontsize(50)
            E0_subplot.yaxis.label.set_fontsize(40)
            n_subplot.xaxis.label.set_fontsize(50)
            n_subplot.yaxis.label.set_fontsize(40)
            Gamma0_subplot.xaxis.label.set_fontsize(50)
            Gamma0_subplot.yaxis.label.set_fontsize(40)
            epsilonB_subplot.xaxis.label.set_fontsize(50)
            epsilonB_subplot.yaxis.label.set_fontsize(40)
            p_subplot.xaxis.label.set_fontsize(50)
            p_subplot.yaxis.label.set_fontsize(40)
            epsilonB_RS_subplot.xaxis.label.set_fontsize(50)
            epsilonB_RS_subplot.yaxis.label.set_fontsize(40)
            p_RS_subplot.xaxis.label.set_fontsize(50)
            p_RS_subplot.yaxis.label.set_fontsize(40)
            theta0_subplot.xaxis.label.set_fontsize(50)
            theta0_subplot.yaxis.label.set_fontsize(40)
            alpha_subplot.xaxis.label.set_fontsize(50)
            alpha_subplot.yaxis.label.set_fontsize(40)
            tprompt_subplot.xaxis.label.set_fontsize(50)
            tprompt_subplot.yaxis.label.set_fontsize(40)
            s_subplot = None
            N_coc_subplot = None
            T_coc_subplot = None
            R_coc_subplot = None
            theta0_coc_subplot = None
            Gamma0_coc_subplot = None

#            return epsilone_subplot, E0_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, alpha_subplot, epsilonB_RS_subplot, p_RS_subplot, epsilone_RS_subplot, tprompt_subplot


    return None,epsilone_subplot,None,epsilone_RS_subplot,None,E0_subplot,n_subplot,s_subplot,Gamma0_subplot,epsilonB_subplot,epsilonB_RS_subplot,p_subplot,None,theta0_subplot,alpha_subplot,None,None,None,None,tprompt_subplot,p_RS_subplot,R_coc_subplot,N_coc_subplot,T_coc_subplot,theta0_coc_subplot,Gamma0_coc_subplot,None,None
