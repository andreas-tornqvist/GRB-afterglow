#This program runs the main fitter

import pymultinest
import numpy as np
import time 
import random
import sys
import os
import glob
import resource
from cosmocalc import cosmocalc



#from specGen import FluxCalc
from fitterFunction import modelFunc

import options
tidPre = time.gmtime()
from decimal import Decimal
global numRuns,surfRingsOut,tdata,FdataInput,errorbarInput,numberOfEmpties
np.set_printoptions(threshold=np.nan)


def myPrior(cube,ndims,nparam):
    #import numpy as np
    cubeIndex = np.array(np.where(parametrar))
    for i in range(ndims): cube[i] = paramLimits[cubeIndex[0,i],0] +  cube[i] * (paramLimits[cubeIndex[0,i],1]-paramLimits[cubeIndex[0,i],0])
    

def logLikelihood(cube,ndims,nparam,cm_FdataInput=None,cm_tdata=None,cm_errorbarInput=None,cm_numberOfEmpties=None,cm_numberOfPoints=None):
    global lowestChi2, numberOfPoints , numRuns , surfRingsOut , tdata , FdataInput , errorbarInput , numberOfEmpties

    if (runOption=='LC') and createMock: FdataInput, tdata, errorbarInput, numberOfEmpties,numberOfEmpties = cm_FdataInput, cm_tdata, cm_errorbarInput, cm_numberOfEmpties, cm_numberOfEmpties
    
    #print loadInputConstants
    #print surfRingsOut
    if runOption == 'fit': 
        for i in range(np.sum(parametrar)): constants[np.array(np.where(parametrar))[0,i]] = cube[i] #Assigning input value from prior
    else:
        startTime = time.time()
        print "Number of parameters: %d"%np.sum(parametrar)
            
    
    if (runOption == 'LC'): 
        #plt.plot(0,0,'b-')
        #plt.plot(0,0,'g-')
        #plt.legend(['B-band','X-rays'])
        
        lightcurve,startJetBreak,endJetBreak,tempGrid = modelFunc(R,constants,fixedRSFSratio,reverseShock,runOption,useEATS,thermalComp,photosphere,printStats,allowPrint,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,createMock,opticalDepth,plotComponents,daysOrSec,fluxAxis,numberOfPoints,printProcess,surfRingsOut,np.sum(parametrar),FS_only,RS_only,thermal_only,chi2_type)
        endTime = time.time()
        print "Timeuse: %f s"%(endTime-startTime)
        return lightcurve , tempGrid
    elif runOption == 'fit': 
        
        if (numRuns == 0) & (loadInputConstants): surfRingsOut = 20
        if (numRuns > 10) and (((numRuns % 10000)>=0) &( (numRuns % 10000) < 10)) :
            runAgain = True
            if numRuns % 10000 == 0: surfRingsOut = 10    #Restart the counter after 10 000 runs, when the fluctuations have decreased
            while runAgain:
                print surfRingsOut
                chi2 , _,_ = modelFunc(R,constants,fixedRSFSratio,reverseShock,runOption,useEATS,thermalComp,photosphere,printStats,allowPrint,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,createMock,opticalDepth,plotComponents,daysOrSec,fluxAxis,numberOfPoints,printProcess,surfRingsOut,np.sum(parametrar),FS_only,RS_only,thermal_only,chi2_type)
                chi22 , _,_ = modelFunc(R,constants,fixedRSFSratio,reverseShock,runOption,useEATS,thermalComp,photosphere,printStats,allowPrint,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,createMock,opticalDepth,plotComponents,daysOrSec,fluxAxis,numberOfPoints,printProcess,surfRingsOut*10,np.sum(parametrar),FS_only,RS_only,thermal_only,chi2_type)
                if np.isinf(chi2) or np.isinf(chi22) or np.isnan(chi2) or np.isnan(chi22): return float('-inf')
                else:compare = np.abs((chi2/chi22)-1)
                if compare > .1:  #Totally unacceptable; run again!
                    surfRingsOut += 20
#                    print "Bad surfRingsOut! compare = %f New surfRingsOut = %d"%(compare,surfRingsOut)
                    if numRuns < 10: numRuns = 0  #If it is during the initial phase, then reset the counter to make ten fresh tests
                    elif ((numRuns % 10000)>=0 & (numRuns % 10000) < 10) : numRuns -= (numRuns % 10000) #Resets the counter
                elif compare > 0.01: #Acceptable for now, but do better next time!
                    surfRingsOut += 10
#                    print "compare = %f. New surfRingsOut = %d"%(compare,surfRingsOut)
                    runAgain = False
                    if numRuns < 10: numRuns = 0  #If it is during the initial phase, then reset the counter to make ten fresh tests
                    elif ((numRuns % 10000)>=0 & (numRuns % 10000) < 10) : numRuns -= (numRuns % 10000) #Resets the counter
                    chi2 , _,_ = modelFunc(R,constants,fixedRSFSratio,reverseShock,runOption,useEATS,thermalComp,photosphere,printStats,allowPrint,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,createMock,opticalDepth,plotComponents,daysOrSec,fluxAxis,numberOfPoints,printProcess,surfRingsOut,np.sum(parametrar),FS_only,RS_only,thermal_only,chi2_type) #Now the number of rings in EATS integrator is corrected. Estimate chi2 again and move on

                else: runAgain = False
        else:   
            chi2 , _,_ = modelFunc(R,constants,fixedRSFSratio,reverseShock,runOption,useEATS,thermalComp,photosphere,printStats,allowPrint,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,createMock,opticalDepth,plotComponents,daysOrSec,fluxAxis,numberOfPoints,printProcess,surfRingsOut,np.sum(parametrar),FS_only,RS_only,thermal_only,chi2_type)

        if printProcess & allowPrint:
        #if outputCounter > nth:
            try:
                if chi2 < lowestChi2:
                    
                    #print 
                    if True: #Set this to True to allow the programme to print the process in a file while running fitting routine
                        try: 
                            interLog = open('interactive_log.txt','r')
                            inText = interLog.read()
                            interLog.close()
                            lowestReadChi2 = (float((inText.split('\n')[1]).split(':')[1]))
                            if lowestReadChi2 < chi2 : #If the file interactive_log.txt has a lower chi2 than the routine, then don't rewrite file. This may occur when running on many cores
                                lowestChi2 = lowestReadChi2
                                numRuns += 1
                                print numRuns
                                return - chi2 / 2
                        except:
                            os.system('touch interactive_log.txt')
                            inText = ''

                    utText =  "New lowest chi2:\n"
                    paramOut = "%s=%s"%(paramNames[0],constants[0])
                    for j in range(1,len(constants)):
                        paramOut = "%s\n%s=%s"%(paramOut,paramNames[j],constants[j])
                    paramOut = "%s\nsurfaceRings = %d"%(paramOut,surfRingsOut)
                    utParams = open('parameters.txt','w')
                    utParams.write(paramOut)
                    utParams.close()
                    
                    printTime = time.gmtime()

                    printParams = "At %d:%d:%d: \nNew lowest chi2: %s\n"%(printTime[3],printTime[4],printTime[5],chi2)
                    for i in range(ndims): 
                        #Adding exceptions for print-out that prints a variable on a different form than the fitter is using
                        printParams = '%s\n%s = %s %s'%(printParams,paramNamesShort[whereParam[i]], 10**constants[whereParam[i]]*(preferredScale[whereParam[i]]=='log') + constants[whereParam[i]]*((preferredScale[whereParam[i]]=='deg')*180/np.pi + (preferredScale[whereParam[i]]=='lin')) , 'degrees' * (preferredScale[whereParam[i]]=='deg'))
#                        if paramNames[whereParam[i]] == 'theta0': printParams =  "%s\ntheta0 = %s degrees"%(printParams,constants[whereParam[i]]*180/np.pi)
#                        elif paramNames[whereParam[i]] == 'alpha': printParams = '%s\nalpha = %s degrees'%(printParams,constants[whereParam[i]]*180/np.pi)
#                        elif paramNames[whereParam[i]] == 'log10(E0)': 
#                            printParams =  "%s\nCollimated energy = %s"%(printParams,(10**constants[whereParam[i]]))
#                        elif paramNames[whereParam[i]] == 'log10(n)': 
#                            printParams = "%s\nn = %s"%(printParams,10**constants[whereParam[i]])
#                        elif paramNames[whereParam[i]] == 't0': 
#                            printParams = "%s\nt0 = %s logt0 = %s"%(printParams,t0,logt0)
#                        elif paramNames[whereParam[i]] == 't_prompt': 
#                            printParams = "%s\nt_prompt = %s s"%(printParams,10**constants[whereParam[i]])
#                        elif paramNames[whereParam[i]] == 'log10(epsilone_FS)':
#                            printParams = "%s\nepsilone_FS = %s"%(printParams,10**constants[whereParam[i]])
#                        elif paramNames[whereParam[i]] == 'log10(epsilone_RS)':
#                            printParams = "%s\nepsilone_RS = %s"%(printParams,10**constants[whereParam[i]])
##                        elif paramNames[whereParam[i]] == 'log10(epsilonB_FS)':
#                            printParams = "%s\nepsilonB_FS = %s"%(printParams,10**constants[whereParam[i]])
#                        elif paramNames[whereParam[i]] == 'log10(epsilonB_RS)':
#                            printParams = "%s\nepsilonB_RS = %s"%(printParams,10**constants[whereParam[i]])
#                        elif paramNames[whereParam[i]] == 'Gamma0':
#                            printParams = "%s\nGamma0 = %s"%(printParams,10**constants[whereParam[i]])
#                        elif paramNames[whereParam[i]] == 'Gamma0 (cocoon)':
#                            printParams = "%s\nGamma0(cocoon) = %s"%(printParams,10**constants[whereParam[i]])

#                        else:
#                            printParams = "%s\n%s = %s"%(printParams,paramNames[whereParam[i]],constants[whereParam[i]])
                    printParams = '%s\nnumRuns = %d'%(printParams,numRuns)
                    printParams = '%s\nsurfaceRings = %d'%(printParams,surfRingsOut)
                    printParams = '%s\nlog10(Likelihood) = %s\n\n'%(printParams,-chi2/2)
                    #print "Field of view started crossing the rim at tobs = %f days and covered the entire rim at tobs = %f days."%(tobs[np.argmin(np.abs(theta-alpha-fovAngle))]/86400.,tobs[np.argmin(np.abs(theta+alpha-fovAngle))]/86400.)
                    lowestChi2 = chi2
                    #print printParams
                    if True: #Set this to True to allow the programme to print the process in a file while running fitting routine
                        interLog = open('interactive_log.txt','w')
                        interLog.write("%s%s"%(printParams,inText))
                        interLog.close()
                        doOnce = False
            except: lowestChi2 = chi2
        numRuns += 1
        return - chi2 / 2


def factor_of_ten(in_value):
    import numpy as np
    return int(np.floor(np.log10(in_value)))
    

def round_off(input_value,round_factor):
    fot = int(np.log10(input_value))
    fot -= (fot < 0)

    return round(input_value / 10**fot,round_factor) * 10**fot

def binner(in_array,n_points,option):
    import numpy as np
    n_bins = len(in_array) / n_points  #Number of bins
    output = np.zeros(n_bins)
    for i in range(n_bins):
        if option == 'linear':
            output[i] = np.sum(in_array[i*n_points:(i+1)*n_points]) / n_points
        elif option == 'square':
            output[i] = np.sqrt(np.sum(in_array[i*n_points:(i+1)*n_points]**2)) / n_points
        else:
            raise NameError('Invalid bin option %s'%option)
    remaining_points = (len(in_array) % n_points)
    if remaining_points != 0: #Number of points is not devidable by number of bins. We have to construct the last bin separatly
        if option == 'linear':
            output = np.append(output,np.sum(in_array[n_bins*n_points:])/remaining_points)
        elif option == 'square':
            output = np.append(output,np.sqrt(np.sum(in_array[n_bins*n_points:]))/remaining_points)
                               
    return output


############################################
#              Loading data....            #
############################################



def loadData():
    global tobsmax, numberOfPoints
#Loading input data

    if useData:     #If useData == True, data is collected from the mock-data directory, otherwise from the real-data dir.
        os.chdir(inputData)
        if allowPrint: print "Getting data files from %s"%inputData
    else: 
        os.chdir(inputData)
        if allowPrint: print "Getting data files from %s"%inputData

    dataFilesTemp = glob.glob("*.%s"%fileExtention)
    outTemp = np.zeros(len(dataFilesTemp))
    tdataTemp,FdataTemp,errorbarTemp,errorIs = np.array([]),np.array([]),np.array([]),np.array([],dtype=bool)
    maxIndex = 0 
    for giveName in range(len(dataFilesTemp)): outTemp[giveName] = float(dataFilesTemp[giveName].split(breakName)[0])

    ### Chosing what lightcurves to load, if user specified the input option 'band=?' ###
    try:
        if bandLimit[1] == 0:   #No upper limit
            bandLimit[1] = 1e30 #A rediculously high limit to facilitate a 'no upper limit'
        inputNumber = len(dataFilesTemp)
        keep = map(bool,np.ones(inputNumber))  #boolean array 'keep' returns True if file should be loaded
        
        for rInd in range(inputNumber):
            if (outTemp[rInd] >= bandLimit[0]) and (outTemp[rInd] <= bandLimit[1]): keep[rInd] = True
            else: keep[rInd] = False

        out = np.zeros(np.sum(keep))
        dataFiles = ['']*np.sum(keep)
        for rKeep in range(inputNumber):
            if keep[rKeep]:
                out[np.sum(keep[:rKeep+1]) - 1] = outTemp[rKeep]
                dataFiles[np.sum(keep[:rKeep+1]) - 1] = dataFilesTemp[rKeep]

    except:
        out = outTemp
        dataFiles = dataFilesTemp

    
    if allowPrint: print "Getting data from %d lightcurves"%len(dataFiles)

    dataIndeces = np.array([0]*(len(dataFiles)+1)) #dataIndeces: The index after the last entry of the input data (se concatenation below)


    numberOfPoints = 0
    for i in range(len(dataFiles)):
        
        #Determining unit in file name
        if inputFreqUnit == 'GHz': out[i] *= 1e9
        elif inputFreqUnit == 'MHz': out[i] *= 1e6
        elif inputFreqUnit == 'cm': out[i] = c / out[i]
    
        try:
            input_data = np.loadtxt(dataFiles[i])
            tTemp = input_data[:,0]
            FRead = input_data[:,1]
            sigmaTemp = input_data[:,2]
            errorIsTemp = np.ones(len(tTemp),dtype=bool)
            if inputTimeFormat == 'day': tTemp *= 86400
            elif inputTimeFormat == 'hour': tTemp *= 3600
            elif inputTimeFormat == 'min': tTemp *= 60


        except:
            a = open(dataFiles[i])
            inText = a.read().split('\n')
            while inText[-1] == '':
                inText = inText[:-1]
            a.close()
        
            if (len(inText[0].split())!=1): #if the first entry in the .dat file is not the number of data points
                inputLength = 0
                for k in range(len(inText)):
                    if (inText[inputLength] != ''):
                        try:
                            map(float,inText[inputLength].split())   #Is the row a row of data points?
                            inputLength += 1
                        except: break
                isFirstTheNumber=0
            else: #if the first entry in the .dat file IS the number of data points
                inputLength = int(inText[0])
                if inputLength != len(inText[1:]): raw_input("Input data is of other length than is stated in top of file!")  #Sanity check
                isFirstTheNumber = 1
        

        
            tTemp,FRead,sigmaTemp,errorIsTemp = np.zeros(inputLength),np.zeros(inputLength),np.zeros(inputLength),np.zeros(inputLength,dtype=bool) #errorIsTemp is an array to say if the errorbar contains true errorbars (True) or upper limits (False)
        #Converting string line to float array
            for j in range(isFirstTheNumber,len(inText)):
                inRow = map(float,inText[j].split())
                tTemp[j-1],FRead[j-1] = np.array(inRow[:2])
            #Converting input time
                if inputTimeFormat == 'day': tTemp[j-1] *= 86400
                elif inputTimeFormat == 'hour': tTemp[j-1] *= 3600
                elif inputTimeFormat == 'min': tTemp[j-1] *= 60

                try: 
                    sigmaTemp[j-isFirstTheNumber] = np.array(inRow[2]) 
                    errorIsTemp[j-isFirstTheNumber] = True
                except: 
                    sigmaTemp[j-isFirstTheNumber] = np.array(inRow[1])   #If the error bar is missing, sigma = flux
                    errorIsTemp[j-isFirstTheNumber] = False
        
            if inputLength != len(tTemp): raw_input("Input data is of other length than is stated in top of file!")  #Sanity check
        
        ### Sorting and binning input data
        indexOrder = np.argsort(tTemp)
        if bin_data_bool:
            tTemp = binner(tTemp[indexOrder],bin_data,'linear')
            FRead = binner(FRead[indexOrder],bin_data,'linear')
            sigmaTemp = binner(sigmaTemp[indexOrder],bin_data,'square')
            errorIsTemp = binner(errorIsTemp,bin_data,'linear')
            indexOrder = np.arange(len(tTemp))
        
        for indSort in range(len(tTemp)):
            tdataTemp = np.append(tdataTemp,tTemp[indexOrder[indSort]])
            FdataTemp = np.append(FdataTemp,FRead[indexOrder[indSort]])
            errorbarTemp = np.append(errorbarTemp,sigmaTemp[indexOrder[indSort]])
            errorIs = np.append(errorIs,errorIsTemp[indexOrder[indSort]])

        inputLength = len(tTemp)
        numberOfPoints += inputLength
        dataIndeces[i+1] = int(dataIndeces[i] + inputLength)
        if (dataIndeces[i+1]-dataIndeces[i]) > maxIndex: maxIndex = dataIndeces[i+1]-dataIndeces[i]
        


    firstColumn,FdataOut,errorbarOut,errorIsOut = np.zeros([len(dataFiles),maxIndex]),np.zeros([len(dataFiles),maxIndex]),np.zeros([len(dataFiles),maxIndex]),np.zeros([len(dataFiles),maxIndex],dtype=bool)

    os.chdir(currentPath)

    for i in range(len(dataFiles)):
        firstColumn[i] = np.concatenate([tdataTemp[dataIndeces[i]:dataIndeces[i+1]],[-1]*(maxIndex-(dataIndeces[i+1]-dataIndeces[i]))])
        FdataOut[i] = np.concatenate([FdataTemp[dataIndeces[i]:dataIndeces[i+1]],[-1]*(maxIndex-(dataIndeces[i+1]-dataIndeces[i]))])
        errorbarOut[i] = np.concatenate([errorbarTemp[dataIndeces[i]:dataIndeces[i+1]],[-1]*(maxIndex-(dataIndeces[i+1]-dataIndeces[i]))])
        errorIsOut[i] = np.concatenate([errorIs[dataIndeces[i]:dataIndeces[i+1]],[-1]*(maxIndex-(dataIndeces[i+1]-dataIndeces[i]))])

        minFdataThis , maxFdataThis = FdataTemp[dataIndeces[i] + np.argmin(FdataTemp[dataIndeces[i]:dataIndeces[i+1]])] , FdataTemp[dataIndeces[i] + np.argmax(FdataTemp[dataIndeces[i]:dataIndeces[i+1]])]
        try: 
            if minFdataThis < minFdata: minFdata = minFdataThis
        except: minFdata = minFdataThis
        try: 
            if maxFdataThis > maxFdata: maxFdata = maxFdataThis
        except: maxFdata = maxFdataThis
        minTdataThis , maxTdataThis = tdataTemp[dataIndeces[i] + np.argmin(tdataTemp[dataIndeces[i]:dataIndeces[i+1]])] , tdataTemp[dataIndeces[i] + np.argmax(tdataTemp[dataIndeces[i]:dataIndeces[i+1]])]
        try: 
            if minTdataThis < minTdata: minTdata = minTdataThis
        except: minTdata = minTdataThis
        try: 
            if maxTdataThis > maxTdata: maxTdata = maxTdataThis
        except: maxTdata = maxTdataThis

        #if plotInput:
        #    plt.figure('plot_input')
        #    print tdataTemp[dataIndeces[i]:dataIndeces[i+1]]
            
            #plt.subplot(int(len(dataFiles)**(1/2.)),int(float(len(dataFiles))/int(len(dataFiles)**(1/2.)))+1,i)
        #    plt.scatter(tdataTemp[dataIndeces[i]:dataIndeces[i+1]]/86400,FdataTemp[dataIndeces[i]:dataIndeces[i+1]],color='%s'%colourCycle[i%len(colourCycle)],s=1)
        #    plt.errorbar(tdataTemp[dataIndeces[i]:dataIndeces[i+1]]/86400,FdataTemp[dataIndeces[i]:dataIndeces[i+1]],yerr=errorbarTemp[dataIndeces[i]:dataIndeces[i+1]],ecolor='%s'%colourCycle[i%len(colourCycle)],linestyle='None')
#    if plotInput:
#        
#        plt.xlim([minTdata/2/86400,maxTdata*2/86400])
#        plt.ylim([minFdata/2,maxFdata*2])
#        plt.xlabel('Time (days)')
#        plt.ylabel('Flux (Jy)')
#        plt.title('GRB %s'%GRBlabel)
#        plt.loglog()

#        if os.path.isdir("Figures/") == False: os.system("mkdir Figures/")
#        #plt.show()
#        plt.savefig("Figures/%s.%s"%(GRBlabel,figTypes))
#        plt.close('plot_input')
#        if allowPrint: print "Figure with input data saved to Figures/%s.%s"%(inputData.split('/')[-2],figTypes)

    

    return firstColumn,FdataOut,errorbarOut,errorIsOut,out


#######################################
#   Analysing and plotting results    #
#######################################

def plot_mstats(x_grid,y_grid_in,n_params,this_plot,gauss_linetype,gauss_line_width,plot_dim,mstats_text,param_number):
    from matplotlib import mlab

    if plot_dim == 'log' or plot_dim == 'lin': x_grid_in = np.copy(x_grid)
    elif plot_dim == 'deg': x_grid_in = x_grid * 180 / np.pi

    mstats_mu = np.zeros(n_params)
    mstats_sigma = np.zeros(n_params)
#    for i_mstats in range(1,len(mstats_text)+1):

    this_line = map(float,mstats_text[param_number+1].split())


    if param_number + 1 != this_line[0]:# Sanity check. If this fails the layout of the 1-stats.dat file is probably changed                                        
        print 'Something occur while reading chains/1-stats.dat. Now exiting.'
        raise SystemExit(0)

    mstats_mu = this_line[1]
    mstats_sigma = this_line[2]

    if len(x_grid_in) == 0: #No points in this distribution. Why is this?
        print 'Parameter %s without distribution! Skipping this one.'%paramNames[whereParam[param_number]]
    else:
        grid_area = np.sum(np.concatenate([[(x_grid_in[1]-x_grid_in[0])/2],(x_grid_in[2:]-x_grid_in[:-2])/2,[(x_grid_in[-1]-x_grid_in[-2])/2]]) * y_grid_in)
    
        xlims = this_plot.get_xlim()
        if plot_dim == 'log': 
            x_plotarray = np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),100)
            plot_pdf = mlab.normpdf(np.linspace(np.log10(min(x_plotarray)),np.log10(max(x_plotarray)),100),mstats_mu,mstats_sigma)
        elif plot_dim == 'lin':
            x_plotarray = np.linspace(xlims[0],xlims[1],100)
            plot_pdf = mlab.normpdf(x_plotarray,mstats_mu,mstats_sigma)
        elif plot_dim == 'deg': 
            x_plotarray = np.linspace(xlims[0],xlims[1],100)
            plot_pdf = mlab.normpdf(x_plotarray,mstats_mu*180/np.pi,mstats_sigma*180/np.pi)
        this_plot.plot(x_plotarray , plot_pdf * grid_area,gauss_linetype,linewidth=gauss_line_width)


def plot_comres(this_plot,fit_pars,color,scale):
    n_comres = np.shape(fit_pars)[0]

    for i_comres in range(n_comres):
        x_lims = this_plot.get_xlim()[0:2]
        y_lims = this_plot.get_ylim()[0:2]

        if fit_pars[i_comres,0] != 0.:
            out_xlims = np.zeros(2)
            new_xlims = np.zeros(2,dtype=bool)
            if fit_pars[i_comres,1] < x_lims[0]: 
                out_xlims[0] = fit_pars[i_comres,1]
                new_xlims[0] = True
            else: out_xlims[0] = x_lims[0]
            if fit_pars[i_comres,2] > x_lims[1]: 
                out_xlims[1] = fit_pars[i_comres,2]
                new_xlims[1] = True
            else: out_xlims[1] = x_lims[1]


            this_plot.plot([fit_pars[i_comres,0],fit_pars[i_comres,0]] , [0.,y_lims[1]*2],color=color[i_comres],linewidth=5.0)
#            this_plot.axvspan(fit_pars[i_comres,1],fit_pars[i_comres,2],alpha=0.5,color=color[i_comres])
            this_plot.plot([fit_pars[i_comres,1],fit_pars[i_comres,1]] , [0., y_lims[1]*2],linestyle='--',color=color[i_comres],linewidth=5.0)
            this_plot.plot([fit_pars[i_comres,2],fit_pars[i_comres,2]] , [0., y_lims[1]*2],linestyle='--',color=color[i_comres],linewidth=5.0)

            if scale == 'log':
                new_xlow = 10**(np.log10(out_xlims[0]) - np.log10(out_xlims[1]/out_xlims[0])/5*new_xlims[0])
                new_xhigh = 10**(np.log10(out_xlims[1]) + np.log10(out_xlims[1]/out_xlims[0])/5*new_xlims[1])
            else:
                new_xlow = out_xlims[0]-(out_xlims[1]-out_xlims[0])/5*new_xlims[0]
                new_xhigh = out_xlims[1] + (out_xlims[1]-out_xlims[0])/5*new_xlims[1]
            this_plot.set_xlim([new_xlow,new_xhigh])
            this_plot.set_ylim(y_lims)


def stats_subplot(this_plot,x_grid_in,y_grid_in,plotLineIn,plot_multinest_gaussian,plot_dim,last_plot,numberOfPlots,param_counter,mstats_text,gauss_linestyle=None,midLinePos=None,midLineColor=None,plot_lims=[float('-inf'),float('inf')]):

    x_grid_selection = np.where((x_grid_in>plot_lims[0]) & (x_grid_in<plot_lims[1]))
    line_width = 5.0
    gauss_line_width = 5.0
    if plot_dim == 'log':
        this_plot.plot(10**x_grid_in[x_grid_selection],y_grid_in[x_grid_selection],plotLineIn,linewidth=line_width)
        if last_plot:
            ylower,yupper = this_plot.get_ylim()
            this_plot.plot(10**np.array([midLinePos,midLinePos]),[0,yupper*1.1],color=midLineColor,linewidth=line_width)
            this_plot.set_ylim([0,yupper])

        
    elif plot_dim == 'lin':
        this_plot.plot(x_grid_in[x_grid_selection],y_grid_in[x_grid_selection],plotLineIn,linewidth=line_width)
        if last_plot:
            ylower,yupper = this_plot.get_ylim()
            this_plot.plot([midLinePos,midLinePos],[0,yupper*1.1],color=midLineColor,linewidth=line_width)
            this_plot.set_ylim([0,yupper])
            deltaMin = (x_grid_in[x_grid_selection][-1] - x_grid_in[x_grid_selection][0])
            facTen = factor_of_ten(x_grid_in[0])
            this_plot.set_xlim([x_grid_in[x_grid_selection][0]-deltaMin,x_grid_in[x_grid_selection][-1]+deltaMin])

    elif plot_dim == 'deg':
        this_plot.plot(x_grid_in[x_grid_selection]*180/np.pi,y_grid_in[x_grid_selection],plotLineIn,linewidth=line_width)
        if last_plot:
            ylower,yupper = this_plot.get_ylim()
            xlower,xupper = this_plot.get_xlim()
            this_plot.plot([180/np.pi*midLinePos,180/np.pi*midLinePos],[0,yupper*1.1],color=midLineColor,linewidth=line_width)
            this_plot.set_xlim([0,xupper])
            this_plot.set_ylim([0,yupper])
        deltaMin = (x_grid_in[x_grid_selection][-1] - x_grid_in[x_grid_selection][0]) * 180 / np.pi
        facTen = factor_of_ten(x_grid_in[x_grid_selection][0])
        lowest,highest = x_grid_in[x_grid_selection][0] * 180 / np.pi , x_grid_in[x_grid_selection][-1] * 180 / np.pi

    if plot_multinest_gaussian: plot_mstats(x_grid_in[x_grid_selection],y_grid_in[x_grid_selection],n_params,this_plot,gauss_linestyle,gauss_line_width,plot_dim,mstats_text,param_counter)    
    for xtick in this_plot.xaxis.get_major_ticks(): xtick.label.set_fontsize(25)
    for ytick in this_plot.yaxis.get_major_ticks(): ytick.label.set_fontsize(25)





def plotter():
    from set_subplots import set_subplots
    from matplotlib import pyplot as plt
    import matplotlib as mpl
    #Analyse output
    print os.path.abspath('.')
    if n_params == 7:
        subHor = 6
        subVer = 3
    
    pathName = '%s/Figures/%s'%(os.path.abspath('..'),GRBlabel)
    if os.path.isdir('%s/Figures/'%(os.path.abspath('..'))) == False: os.system('mkdir ../Figures/')
    if os.path.isdir('%s/Figures/%s/'%(os.path.abspath('..'),GRBlabel)) == False: os.system('mkdir ../Figures/%s/'%GRBlabel)    
    os.system('mkdir %s'%(pathName))
    mpl.rc('xtick', labelsize=40) 
    mpl.rc('ytick', labelsize=40) 
    #axis_font = {'fontname':'Arial', 'size':'14'}
#    plt.title('Probability plot of %s'%GRBlabel,fontsize=40)
    getLimitArray = np.zeros(n_params)

    #Number of plots
    numberOfPlots = input('Number of plots: ')

    #Create average plots
    sum_plots_input = raw_input('Do you want to create average plots? (y/[n]): ')
    if sum_plots_input == 'y' or sum_plots_input == 'Y': sum_plots = True
    else: sum_plots = False

    compare_results_input = raw_input('Do you want to compare with and plot other fit parameters? ([y]/n): ')
    compare_results = (compare_results_input=='Y') or (compare_results_input=='y') or (compare_results_input=='')

    #Print input parameter line?
    if not compare_results:
        printInputLineInput = raw_input('Plot an input value line (use for fitting mock observations) ([y]/n):')
        printInputLine = (printInputLineInput == '') or (printInputLineInput == 'y') or (printInputLineInput == 'Y')
    else:
        n_comres = input('Number of models to input and compare with: ')
        printInputLine = False
    if printInputLine: midLineColor = raw_input('Color of the input value line: ')
    else: midLineColor = None
    pathToPlot = ['']*numberOfPlots        
    plotLine = np.zeros(numberOfPlots,dtype='S10')
    currentPlotPath = os.path.abspath('.')
    #n_average_plots = np.ones(numberOfPlots,dtype=int)
    if compare_results:
        compare_input = np.zeros([n_comres,3,n_params])
        comres_infile = np.zeros(n_comres,dtype='S256')
        comres_color = np.zeros(n_comres,dtype='S20')

    for plotI in range(numberOfPlots):
        if sum_plots:
            #n_average_plots[plotI] = input('Number of plots to average over in plot %d: '%(plotI+1))
            pathToPlot[plotI] = raw_input('Type general path (using asterix) to the %d%s location of all chains folders: '%(plotI+1,'st'*(plotI==0)+'nd'*(plotI==1)+'rd'*(plotI==2)+'th'*(plotI>2)))
            plot_multinest_gaussian = False
            gauss_color = ''
            gauss_linetype = ''
        else:
            pathToPlot[plotI] = os.path.abspath(raw_input('Path where to find chains/ no. %d: '%(plotI+1)))
            plot_multinest_gaussian_input = raw_input('Plot MultiNest output parameters as gaussian? ([y]/n): ')
            plot_multinest_gaussian = plot_multinest_gaussian_input == 'y' or plot_multinest_gaussian_input == 'Y' or plot_multinest_gaussian_input == ''

            try:gauss_linetype
            except: 
                gauss_linetype = np.zeros(numberOfPlots,dtype='S20')
                gauss_color = np.zeros(numberOfPlots,dtype='S20')
            if plot_multinest_gaussian: 
                gauss_linetype[plotI] = np.copy(raw_input('Gaussian linetype: ([-],--,-.,:): '))
                if gauss_linetype[plotI] == '': gauss_linetype = '-'
                gauss_color[plotI] = np.copy(raw_input('Gaussian color: (default=g) '))
                if gauss_color[plotI] == '': gauss_color[plotI] = 'g'
            
        plotLine[plotI] = raw_input('Line color. Default=\'b\': ')
        if plotLine[plotI] == '': plotLine[plotI] = r'b'
        
    if compare_results:   #Input fit parameters for comparison
        for i_comps in range(n_comres):
            comres_infile[i_comps] = raw_input('Name of file with fit parameters%s. Leave empty for manual input: '%((' for model %d'%(plotI+1))*(numberOfPlots>1)))
            comres_color[i_comps] = raw_input('Color: ')
            if comres_infile[i_comps] == '':   #Manual input
                for i_comres in range(n_params):
                    print '\nEnter fit parameters of model %d. Leave empty to skip the parameter. %s'%(i_comps+1,'Enter 0 as 1-sigma to not plot the uncertainty range')
                    try:compare_input[i_comps,0,i_comres] = input('Central value for %s: '%(paramNames[whereParam[i_comres]]))
                    except: 
                        compare_input[i_comps,0,i_comres] = 0.
                        compare_input[i_comps,1,i_comres] = 0.
                        compare_input[i_comps,2,i_comres] = 0.
                    else:
                        compare_input[i_comps,1,i_comres] = input('Lower 1-sigma for %s: '%(paramNames[whereParam[i_comres]]))
                        if compare_input[i_comps,1,i_comres] != 0:  #Will not plot uncertainty range if <-- == 0
                            compare_input[i_comps,2,i_comres] = input('Upper 1-sigma for %s: '%(paramNames[whereParam[i_comres]]))
                        else:
                            compare_input[i_comps,1,i_comres] = compare_input[i_comps,0,i_comres]
                            compare_input[i_comps,2,i_comres] = compare_input[i_comps,0,i_comres]
                        if compare_input[i_comps,1,i_comres] < 0: #Transforms relative standard dev. to absolute
                            compare_input[i_comps,1,i_comres] += compare_input[i_comps,0,i_comres]
                            compare_input[i_comps,2,i_comres] += compare_input[i_comps,0,i_comres]

                save_comres_filename = raw_input('File name to save fit parameters of model %d: '%(i_comps+1))
                np.savetxt(save_comres_filename,compare_input[i_comps])

            else:
                try:
                    comres_load_file = np.loadtxt(comres_infile[i_comps])
                    compare_input[i_comps] = np.copy(comres_load_file)
                except:
                    if not os.path.isfile(comres_infile[i_comps]): print 'File %s does not exist'%comres_infile[i_comps]
                    elif np.shape(np.loadtxt(comres_infile[i_comps])) != np.shape(compare_input[i_comps]): print 'Number of parameters in file %s is not consistent with number parameters specified in options.py'
                    print 'Now exiting.'
                    raise SystemExit(0)
    
            
            

    linePlotScale = ['']*n_params
    create_temp_catalogue = False

    #Setting subplots

    list_subplots = set_subplots(reverseShock,n_params,parametrar)#[None,epsilone_subplot,None,epsilone_RS_subplot,None,E0_subplot,n_subplot,s_subplot,Gamma0_subplot,epsilonB_subplot,epsilonB_RS_subplot,p_subplot,None,theta0_subplot,alpha_subplot,None,None,None,None,tprompt_subplot,p_RS_subplot,R_coc_subplot,N_coc_subplot,T_coc_subplot,theta0_coc_subplot,Gamma0_coc_subplot,None,None]
    print len(list_subplots)
    print len(parametrar)

#    if reverseShock:
#        if n_params == 10: epsilone_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, alpha_subplot, epsilonB_RS_subplot, p_RS_subplot,# epsilone_RS_subplot, tprompt_subplot = set_subplots(reverseShock,n_params)
#        elif n_params == 11: epsilone_subplot, E0_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, epsilonB_RS_subplot, p_RS_subplot, epsilone_RS_subplot, tprompt_subplot = set_subplots(reverseShock,n_params)
#        elif n_params == 12: epsilone_subplot, E0_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, alpha_subplot, epsilonB_RS_subplot, p_RS_subplot, epsilone_RS_subplot, tprompt_subplot = set_subplots(reverseShock,n_params)

#    else:
#        if thermalComp:
#            epsilone_subplot, E0_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, alpha_subplot , T_coc_subplot , R_coc_subplot , Gamma0_subplot = set_subplots(reverseShock,n_params)
#        epsilone_subplot, E0_subplot, n_subplot, Gamma0_subplot, epsilonB_subplot, p_subplot, theta0_subplot, alpha_subplot = set_subplots(reverseShock,n_params)

#    list_subplots[whereParam] = 



    lowerLims,upperLims = np.zeros(n_params) , np.zeros(n_params)

    for j in range(numberOfPlots):
        if sum_plots: sum_average_paths = glob.glob(pathToPlot[j])
        if sum_plots and (len(sum_average_paths) > 1): 
            
            x_grid_lims = np.zeros([n_params,2])
            for sum_average in range(len(sum_average_paths)):
                os.chdir(currentPlotPath)
                os.chdir(sum_average_paths[sum_average])
        
                a = pymultinest.Analyzer(n_params = n_params)
                s = a.get_stats()
                p = pymultinest.PlotMarginalModes(a)
                i_counter = 0
                if os.path.isdir('%s/temp_catalogue'%currentPlotPath) and j == 0 and sum_average == 0: 
                    os.system('rm %s/temp_catalogue/*'%currentPlotPath)
                elif not os.path.isdir('%s/temp_catalogue'%currentPlotPath): 
                    create_temp_catalogue = True
                    os.system('mkdir %s/temp_catalogue'%currentPlotPath)
                for i in np.where(parametrar)[0]:

                    x_grid,y_grid = p.plot_marginal(i_counter, with_ellipses = True, with_points = False, grid_points=50)
                    np.savetxt('%s/temp_catalogue/grid_%d%d.txt'%(currentPlotPath,i,sum_average),[x_grid,y_grid])

                    #Finding extremes for the x-axis
                    if sum_average == 0:
                        x_grid_lims[i_counter] = np.array([min(x_grid),max(x_grid)])
                    else:
                        if x_grid_lims[i_counter,0] > min(x_grid): x_grid_lims[i_counter,0] = np.copy(min(x_grid))
                        if x_grid_lims[i_counter,1] < max(x_grid): x_grid_lims[i_counter,1] = np.copy(max(x_grid))
                    i_counter += 1
            os.chdir(currentPlotPath)
            mstats_text = None
        else: 
            
                                                            

            os.chdir(pathToPlot[j])
        
            a = pymultinest.Analyzer(n_params = n_params)
            s = a.get_stats()
            p = pymultinest.PlotMarginalModes(a)

            #Reading mutlinest output from chains/1-stats.dat
            if plot_multinest_gaussian:
                open_mstats = open('chains/1-stats.dat')
                read_mstats = open_mstats.read().split('Dim No.')
                open_mstats.close()
                mstats_text = read_mstats[1].split('\n')
            else:
                mstats_text = 'None'
        
        limCounter = 0


        if True:
            for i in np.where(parametrar)[0]:
            ### Getting data form modified pymultinest.plot module

                if sum_plots and (len(sum_average_paths) > 1): #Creating average x-and-y grid
                    n_ave = 100
                    x_grid = np.linspace(x_grid_lims[limCounter,0],x_grid_lims[limCounter,1],n_ave)
                    y_grid = np.zeros(n_ave)
                    y_counter = 0
                    for average_y in range(len(sum_average_paths)):
                        read_x_grid, read_y_grid = np.loadtxt('%s/temp_catalogue/grid_%d%d.txt'%(currentPlotPath,i,average_y))
                        for y_index in range(n_ave):
                            #Interpolating
                            nearest_under = np.argmin(np.abs(read_x_grid - x_grid[y_index])) #Finding the index in the read-in x_grid closest to the active x_grid element
                            nearest_under -= (read_x_grid[nearest_under] > x_grid[y_index]) #Making sure we have the index just under the active x_grid index
                            if nearest_under < 0: continue
                            elif nearest_under > (len(read_y_grid)-2): break
                            y_grid[y_index] += (read_y_grid[nearest_under] * (read_x_grid[nearest_under+1] - x_grid[y_index]) + read_y_grid[nearest_under+1] * (x_grid[y_index] - read_x_grid[nearest_under])) / (read_x_grid[nearest_under+1] - read_x_grid[nearest_under])  #Interpolation to find y_grid value on the new grid
                            
                        y_counter += 1
                    y_grid /= y_counter   #Averaging

                else:
                    x_grid,y_grid = p.plot_marginal(limCounter, with_ellipses = True, with_points = False, grid_points=50)



            #x_grid[x_grid<0.01] = 0.
                n_minlim = np.max(y_grid) * 0.01
                lowest_n , highest_n = 0,len(x_grid)-1
                while y_grid[lowest_n] < n_minlim: lowest_n += 1
                while y_grid[highest_n] < n_minlim: highest_n -= 1
                try:
                    if lowest_n != 0:x_grid,y_grid = x_grid[lowest_n-1:highest_n+2] , y_grid[lowest_n-1:highest_n+2]
                    else: x_grid,y_grid = x_grid[0:highest_n+2] , y_grid[0:highest_n+2]
                except: #In case highest_n is the last index
                    if lowest_n != 0:x_grid,y_grid = x_grid[lowest_n-1:highest_n+1] , y_grid[lowest_n-1:highest_n+1]
                    else: x_grid,y_grid = x_grid[0:highest_n+1] , y_grid[0:highest_n+1]

#                if not reverseShock:  #Forward shock (no reverse shock)
#                if i == 1: #epsilon_e
#                    this_plot = epsilone_subplot

#                elif i == 3 and reverseShock:
#                    this_plot = epsilone_RS_subplot
#                    stats_subplot(epsilone_RS_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'log',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor,[float('-inf'),0])
#                elif i == 5: #E0
#                    this_plot = E0_subplot
#                    stats_subplot(E0_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'log',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor)
#                elif i == 6: #n
#                    this_plot = n_subplot
#                    stats_subplot(n_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'log',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor)
#                elif i == 8: #Gamma0
#                    this_plot = Gamma0_subplot
#                    stats_subplot(Gamma0_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'log',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor)
#                elif i == 9:#epsilon_B
#                    this_plot = epsilonB_subplot
#                    stats_subplot(epsilonB_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'log',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor,[float('-inf'),0])
#                elif i == 10 and reverseShock:
#                    this_plot = epsilonB_RS_subplot
#                    stats_subplot(epsilonB_RS_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'log',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor,[float('-inf'),0])
#                elif i == 11: #p
#                    this_plot = p_subplot
#                    stats_subplot(p_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'lin',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor)
#                elif i == 13:# theta0
#                    this_plot = theta0_subplot
#                    stats_subplot(theta0_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'degree',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor,[0.,float('inf')])
#                elif i == 14:# alpha
#                    this_plot = alpha_subplot
#                    stats_subplot(alpha_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'degree',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor,[0.,float('inf')])
#                elif i == 19 and reverseShock: #t_prompt
#                    this_plot = tprompt_subplot
#                    stats_subplot(tprompt_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'log',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor,[0.,float('inf')])
#                elif i == 20 and reverseShock: #p_RS
#                    this_plot = p_RS_subplot
#                    stats_subplot(p_RS_subplot,x_grid,y_grid,plotLine[j],plot_multinest_gaussian,'%s%s'%(gauss_color[j],gauss_linetype[j]),'lin',(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,constants[i],midLineColor)
                if plot_multinest_gaussian: stats_subplot(list_subplots[i],x_grid,y_grid,plotLine[j],plot_multinest_gaussian,preferredScale[i],(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,'%s%s'%(gauss_color[j],gauss_linetype[j]),constants[i],midLineColor,paramLimits[i])                    
                else: stats_subplot(list_subplots[i],x_grid,y_grid,plotLine[j],plot_multinest_gaussian,preferredScale[i],(j==(numberOfPlots-1)) and printInputLine,numberOfPlots,limCounter,mstats_text,None,constants[i],midLineColor,paramLimits[i])                    

                
                if compare_results:   #Plot input fit parameters
                    plot_comres(this_plot,compare_input[:,:,limCounter],comres_color,preferredScale[i])
    
                limCounter += 1
    #Deleting temp_catalogue
    if create_temp_catalogue: os.system('rm %s/temp_catalogue -r'%currentPlotPath)
        ### Setting ticsk

        
        
    
        #plt.ylim([0,getLimitArray[i]*1.1])
        #print "ylim set"
    os.chdir(currentPlotPath)

    
    plt.tight_layout()
    if os.path.isfile('%s/probability_%s.%s'%(pathName,GRBlabel,figTypes)):
        del_fig_file = raw_input('Overwrite probility_%s.%s? ([y]/n): '%(GRBlabel,figTypes))
        if (del_fig_file == 'n') or (del_fig_file == 'N'): 
            output_filename = raw_input('Type file name without suffix: ')
        else:output_filename = 'probability_%s'%(GRBlabel)
    else:output_filename = 'probability_%s'%(GRBlabel)
    plt.savefig('%s/%s.%s'%(pathName,output_filename,figTypes))
    plt.close()
    #os.system('pdftk %s/*.pdf cat output %s/merged.pdf'%(pathName,pathName))
    if allowPrint: print "Probablility plots saved in %s/%s.%s"%(pathName,output_filename,figTypes)



### This function gathers the information from the chains/1-stats.dat file and prints an analysis of them compared to the input value
def gather_multinest_stats():
    stats_path_input = raw_input('Path where to find the chains folder containing the 1-stats.dat file. Type general path using asterix, separated by space to analyse many or a file containing a list of all paths separated by newlines: ')
    write_tex_files_input = raw_input('Write tex-files? ([y]/n): ')
    write_tex_files = (write_tex_files_input == '') or (write_tex_files_input == 'y') or (write_tex_files_input == 'Y')
    if not write_tex_files:
        read_output_matrix_input = raw_input('Analyse written files? ([y]/n): ')
        read_output_matrix = (read_output_matrix_input == '') or (read_output_matrix_input == 'Y') or (read_output_matrix_input == 'y')
    
    if len(stats_path_input.split()) > 1: ### Many inputs
        stats_path = np.zeros(0,dtype='S256')
        for i in range(len(stats_path_input.split())):
            stats_path = np.append(stats_path,glob.glob('%s/chains/1-stats.dat'%stats_path_input.split()[i]))
    else:
       try: 
           open_stats_path = open(stats_path_input,'r')
           stats_path = open_stats_path.read().split('\n')
           open_stats_path.close()
           while stats_path[-1] == '': stats_path = stats_path[:-1]
       except:
           stats_path = glob.glob('%s/chains/1-stats.dat'%stats_path_input)
    stats_path_to_chains = np.zeros(np.shape(stats_path),dtype='S256')

    true_input = constants[whereParam]    

    lower_limit , upper_limit = np.zeros([len(stats_path),n_params],dtype=bool) , np.zeros([len(stats_path),n_params],dtype=bool)

    dist_type = np.zeros((len(stats_path),n_params),dtype='S40')
    type_of_distribution = np.zeros((len(stats_path),n_params),dtype=int)
    
    stats_gaussian = np.zeros([len(stats_path),n_params,2])
    stats_max_like = np.zeros([len(stats_path),n_params])
#    stats_MAP = np.zeros([len(stats_path),n_params,2])
        
    sigma_percentage = np.zeros([len(stats_path),n_params])
    sigma_percentage_text = np.zeros([len(stats_path),n_params],dtype='S20')
    true_input_array = np.ones([len(stats_path),n_params]) * true_input
    mu_diff = np.zeros([len(stats_path),n_params])
    mu_diff_sign = np.zeros([len(stats_path),n_params],dtype=bool)

    bestfit_diff = np.zeros([len(stats_path),n_params])
    bestfit_diff_sign = np.zeros([len(stats_path),n_params],dtype=bool)
    same_side_sign = np.zeros([len(stats_path),n_params],dtype=bool)

    diff_sign_text = np.zeros(np.shape(same_side_sign),dtype='S20')
    inside_sigma = np.zeros([len(stats_path),n_params],dtype=bool)
    inside_sigma_text = np.zeros(np.shape(inside_sigma),dtype='S20')

#    write_cover_file_input = raw_input('Create cover file? ([y]/n): ')
#    write_cover_file = write_cover_file_input == '' or write_cover_file_input == 'y' or write_cover_file_input == 'Y'
#    if write_cover_file: cover_text = ''

    #if read_output_matrix:
        ### Creating a 3-D tensor with dimensions Number of files x Number of entries x Number of parameters
    #    input_data = np.zeros([1,6,n_params])

    for i_path in range(len(stats_path)):
        if write_tex_files:
            outtext = '\\begin{tabular}{l|c|c|c|c|c|c|c|c}\nParameter & Dist. type & T.I. loc. & Max. Like. & $\mu$ & T.I. & No. S.Ds. & Bias & Size of $\sigma$\\\\\n'
        stats_path_to_chains[i_path] = os.path.abspath('%s'%('/'.join(stats_path[i_path].split('/')[:-2])))

        if write_tex_files: os.chdir(stats_path_to_chains[i_path])

        path_name_split = os.path.abspath(stats_path_to_chains[i_path]).split('/')
        if len(path_name_split[-1]) > 4: path_name_split[-1] = path_name_split[-1][:4]


        if write_tex_files:
            output_filename = raw_input('Output file name without suffix [%s_%s_%s.tex]: '%(path_name_split[-3],path_name_split[-2],path_name_split[-1]))
        else:
            output_filename = '%s_%s_%s'%(path_name_split[-3],path_name_split[-2],path_name_split[-1])
        if output_filename == '': output_filename = '%s_%s_%s'%(path_name_split[-3],path_name_split[-2],path_name_split[-1])
        if os.path.isfile('%s.tex'%output_filename) and write_tex_files:
            ask_del_input = raw_input('Delete existing file %s.tex? ([y]/n): '%output_filename)
            ask_del = ask_del_input == 'y' or ask_del_input == 'Y' or ask_del_input == ''
            if not ask_del:
                print 'Will not overwrite file %s.tex. Now exiting.'%output_filename
                raise SystemExit(0)
            else: os.system('rm  %s.tex'%output_filename)


        ####################################
        ### Loading written matrix files ###
        ####################################

        if read_output_matrix:
            input_matrix_files = glob.glob('%s/%s_mode*.dat'%(stats_path_to_chains[i_path].split('chains/1-stats.dat')[0] , output_filename))

            for open_file_modes in range(len(input_matrix_files)):  ### Looping over all modes
                try:
                    input_data_temp = np.copy(input_data)
                    input_data = np.zeros([np.shape(input_data)[0]+1,np.shape(input_data)[1],np.shape(input_data)[2]])
                    input_data[-1] = np.loadtxt(input_matrix_files[open_file_modes])
                    input_data_path = np.append(input_data_path,'%s'%input_matrix_files[open_file_modes])
                    input_data_model_list = np.append(input_data_model_list,np.array(path_name_split[-1]))
                except: 
                    input_data = np.zeros([1,6,n_params])
                    input_data_path = np.zeros(1,dtype='S256')
                    input_data_path[0] = '%s'%input_matrix_files[open_file_modes]
                    input_data_model_list = np.array([path_name_split[-1]])
                    try:input_data[0] = np.loadtxt(input_matrix_files[open_file_modes])
                    except:
                        print 'Wrong number of parameters specified in options.py. Now exiting.'
                        raise SystemExit(0)
            
            

            continue

        ######################################
        ### Loading posterior distribution ###
        ######################################

        
        try:
            a = pymultinest.Analyzer(n_params = n_params)
            s = a.get_stats()
            p = pymultinest.PlotMarginalModes(a)
        except:
            print 'Could not read files in %s. Now exiting'%stats_path_to_chains[i_path]
            raise SystemExit(0)

        for read_x_grid in range(n_params):

            x_grid,y_grid = p.plot_marginal(read_x_grid, with_ellipses = True, with_points = False, grid_points=50)
            x_grid_selection = np.where((x_grid>paramLimits[whereParam[read_x_grid],0]) & (x_grid<paramLimits[whereParam[read_x_grid],1])) #Reducing grid to only contain points within user specified prior range
            x_grid,y_grid = x_grid[x_grid_selection] , y_grid[x_grid_selection]


            ### Detecting lower and upper limits ;  type_of_distribution: see note from 12/8 -15
            
            if y_grid[0] > (max(y_grid)*0.5): 
                upper_limit[i_path,read_x_grid] = True
                dist_type[i_path,read_x_grid] = 'UL'
                type_of_distribution[i_path,read_x_grid] = 2
            elif y_grid[0] > (max(y_grid)*0.2):
                upper_limit[i_path,read_x_grid] = True
                dist_type[i_path,read_x_grid] = 'PUL'
                type_of_distribution[i_path,read_x_grid] = 3
            if y_grid[-1] > (max(y_grid)*0.5): 
                lower_limit[i_path,read_x_grid] = True
                dist_type[i_path,read_x_grid] = 'LL'
                type_of_distribution[i_path,read_x_grid] = 4
            elif y_grid[-1] > (max(y_grid)*0.2): 
                lower_limit[i_path,read_x_grid] = True
                dist_type[i_path,read_x_grid] = 'PLL'
                type_of_distribution[i_path,read_x_grid] = 5
            if upper_limit[i_path,read_x_grid] == False and lower_limit[i_path,read_x_grid] == False:                 
                dist_type[i_path,read_x_grid] = '{\color{green} CD }'
                type_of_distribution[i_path,read_x_grid] = 0
            if upper_limit[i_path,read_x_grid] and lower_limit[i_path,read_x_grid]: 
                dist_type[i_path,read_x_grid] = '{\\color{red} UC}'
                type_of_distribution[i_path,read_x_grid] = 1
                





        ### Loading multinest outputs from chains/1-stats.dat
        os.chdir(currentPath)
        open_stats = open(os.path.abspath(stats_path[i_path]))
        read_stats = open_stats.read()
        open_stats.close()

        number_of_modes = int(read_stats.split('Total Modes Found:')[1].split()[0])
        print '\n%s has %d mode%s'%(os.path.abspath(stats_path[i_path]),number_of_modes,'s'*(number_of_modes>1))

        

        ### Reading in the stats ###
        local_log_evidence = np.zeros([number_of_modes,n_params])
        global_log_evidence = np.zeros(n_params)

        global_log_evidence[0] = float(read_stats.split('+/-')[0].split()[-1])
        global_log_evidence[1] = float(read_stats.split('+/-')[1].split()[0])

        #if number_of_modes == 1:
        for i_modes in range(number_of_modes):
            if number_of_modes > 1 and write_tex_files:
                outtext += 'mode %d &&&&&\\\\\n'%(i_modes+1)
            stats_split = read_stats.split('Dim No.')
            local_log_evidence[i_modes,0] = float(read_stats.split('Local Log-Evidence')[i_modes+1].split()[0])
            local_log_evidence[i_modes,1] = float(read_stats.split('Local Log-Evidence')[i_modes+1].split()[2])
            
            

                                                                                                           
                                                
#            print stats_split
            for i_stats in range(n_params):
                
                
                stats_gaussian[i_path,i_stats] = map(float,stats_split[3*i_modes+1].split('\n')[i_stats+1].split())[1:3]
                
#                print map(float,stats_split[1].split('\n')[i_stats+1].split())[1:3]
                
                stats_max_like[i_path,i_stats] = map(float,stats_split[3*i_modes+2].split('\n')[i_stats+1].split())[1]
#                stats_MAP[i_path,i_stats] = map(float,stats_split[number_of_modes*i_modes+3].split('\n')[i_stats+1].split())
#            if np.sum(stats_max_like==stats_MAP) == np.size(stats_max_like):
#                print ''#\nMaximum Likelihood Parameters equal the MAP Parameters'
#            else:
#                print ''#\n\n\nMaximum Likelihood Paramters does not equal the MAP Parameters! What to do now?\n\n\n'
                
                ### Percentage of the sigma/mu 
                ### ### For log-scale, print how many decades
                ### ### For lin-scale, print the percentage
        
                if preferredScale[whereParam[i_stats]] == 'lin' or preferredScale[whereParam[i_stats]] == 'deg': #Parameter is fitted on its linear
                    if preferredScale[whereParam[i_stats]] == 'deg': value_factor = 180/np.pi  #Factor to write out degrees in output table
                    else: value_factor = 1.
                    maxlike_exp = int(np.log10(stats_max_like[i_path,i_stats]*value_factor)) * ((np.log10(stats_max_like[i_path,i_stats]*value_factor) > 3) or (np.log10(stats_max_like[i_path,i_stats]*value_factor) < -3))
                    maxlike_base = 10**(np.log10(stats_max_like[i_path,i_stats]*value_factor) - maxlike_exp)  #Gives power-of-ten base if too large or too small to write out
                    
                    mu_exp = int(np.log10(stats_gaussian[i_path,i_stats,0]*value_factor)) * ((np.log10(stats_gaussian[i_path,i_stats,0]*value_factor) > 3) or (np.log10(stats_gaussian[i_path,i_stats,0]*value_factor) < -3))
                    mu_base = 10**(np.log10(stats_gaussian[i_path,i_stats,0]*value_factor)-mu_exp)  #Gives power-of-ten base if too large or too small to write out

                    trueinput_exp = int(np.log10(true_input[i_stats]*value_factor)) * ((np.log10(true_input[i_stats]*value_factor) > 3) or (np.log10(true_input[i_stats]*value_factor) < -3))
                    trueinput_base = 10**(np.log10(true_input[i_stats]*value_factor) - trueinput_exp)


                    if stats_gaussian[i_path,i_stats,1] == 0.: sigma_percentage_text[i_path,i_stats] = 'inf'
                    else:
                        sigma_percentage[i_path,i_stats] = stats_gaussian[i_path,i_stats,1] / stats_gaussian[i_path,i_stats,0] * 100
                        sigma_percentage_text[i_path,i_stats] = '%s \\%% '%(round_off(sigma_percentage[i_path,i_stats]*100,2))
#                elif preferredScale[whereParam[i_stats]] == 'deg':
#                    if stats_gaussian[i_path,i_stats,1] == 0.: sigma_percentage_text[i_path,i_stats] = 'inf'
#                    else:
#                        sigma_percentage_text[i_path,i_stats] = '%s $^{\circ}$'%(round_off( stats_gaussian[i_path,i_stats,1]*180/np.pi ,2))
                elif preferredScale[whereParam[i_stats]] == 'log':
                    maxlike_exp = int(stats_max_like[i_path,i_stats]) * ((stats_max_like[i_path,i_stats] > 3) or (stats_max_like[i_path,i_stats] < -3))
                    maxlike_base = 10**(stats_max_like[i_path,i_stats] - maxlike_exp)  #Gives power-of-ten base if too large or too small to write out
                    
                    mu_exp = int(stats_gaussian[i_path,i_stats,0]) * ((stats_gaussian[i_path,i_stats,0] > 3) or (stats_gaussian[i_path,i_stats,0] < -3))
                    mu_base = 10**(stats_gaussian[i_path,i_stats,0]-mu_exp)  #Gives power-of-ten base if too large or too small to write out

                    trueinput_exp = int(true_input[i_stats]) * ((true_input[i_stats] > 3) or (true_input[i_stats] < -3))
                    trueinput_base = 10**(true_input[i_stats] - trueinput_exp)

                    if stats_gaussian[i_path,i_stats,1] == 0.: sigma_percentage_text[i_path,i_stats] = 'inf'
                    else:
                        sigma_percentage[i_path,i_stats] = stats_gaussian[i_path,i_stats,1]   #How many decades does the standard deviation cover?
                        sigma_percentage_text[i_path,i_stats] = '%s dec'%(((-1)**(sigma_percentage[i_path,i_stats]<0))*round_off(np.abs(sigma_percentage[i_path,i_stats]),2))
                else:
                    print 'Scale is not defined for parameter %s. Now exiting.'%(paramNames[whereParam[i_stats]])
                    raise SystemExit(0)




#        else:
#            raw_input('\n\n\nMore than one mode nb nb nb!!!! Now write the code to analyse this file!\n\n\n')
                             
    ####################################################
    ### Comparing multinest output to the true input ###
    ####################################################

    

    ### Difference between true input and mean value. Negative - over estimated; positive - under estimated
                mu_diff[i_path,i_stats] = true_input[i_stats] - stats_gaussian[i_path,i_stats,0]
                mu_diff_sign[i_path,i_stats] = mu_diff[i_path,i_stats] > 0
    ### Difference between true input and the best-fit value. Negative - over estimated; positive - under estimated
                bestfit_diff[i_path,i_stats] = true_input[i_stats] - stats_max_like[i_path,i_stats]
                bestfit_diff_sign[i_path,i_stats] = bestfit_diff[i_path,i_stats] > 0

                same_side_sign[i_path,i_stats] = mu_diff_sign[i_path,i_stats] == bestfit_diff_sign[i_path,i_stats] #Are mu and the best-fit on the same side of the true input?

    

#                where_diff_sign[i_path,i_stats] = np.where(same_side_sign[i_path,i_stats])
#                where_diff_sign_false = np.where(same_side_sign==False)

                if same_side_sign[i_path,i_stats]:
#                for i_bestfit in range(len(where_diff_sign[0])):
                    diff_sign_text[i_path,i_stats] = mu_diff_sign[i_path,i_stats] * '{\\color{red}UE}' + (mu_diff_sign[i_path,i_stats]==False) * '{\\color{blue}OE}'
                else:
# i_bestfit_false in range(len(where_diff_sign_false[0])):

                    diff_sign_text[i_path,i_stats] = '$\mu$%sb%s'%(mu_diff_sign[i_path,i_stats] * 'U' + (mu_diff_sign[i_path,i_stats]==False) * 'O'  ,  bestfit_diff_sign[i_path,i_stats] * 'U' + (bestfit_diff_sign[i_path,i_stats]==False) * 'O')
    
    ### Is the true input inside the standard deviation?
                inside_sigma[i_path,i_stats] = np.array((true_input[i_stats] > (stats_gaussian[i_path,i_stats,0] - stats_gaussian[i_path,i_stats,1])) & (true_input[i_stats] < (stats_gaussian[i_path,i_stats,0] + stats_gaussian[i_path,i_stats,1])))


                if inside_sigma[i_path,i_stats]: inside_sigma_text[i_path,i_stats] = '{\color{green} IS}'
                else: inside_sigma_text[i_path,i_stats] = '{\color{red} OS}'

 
   


    
    #######################
    ### Printing output ###
    #######################


                no_s_d_ =  round_off(np.abs(mu_diff[i_path,i_stats]),2)
                
                maxlike_text = '$%s%s%s$'%(round_off(maxlike_base,2) , (' \\times 10^{%d}'%maxlike_exp)*(maxlike_exp != 0) , ' ^{\circ}'*(preferredScale[whereParam[i_stats]]=='deg')) #Only writes out the power of ten if it is needed
                
                mu_text = '$%s%s%s$'%(round_off(mu_base,2) , (' \\times 10^{%d}'%mu_exp)*(mu_exp!=0) , ' ^{\circ}'*(preferredScale[whereParam[i_stats]]=='deg')) #Only writes out the power of ten if it is needed
                
                trueinput_text = '$%s%s%s$'%(round_off(trueinput_base,2) , (' \\times 10^{%d}'%trueinput_exp)*(trueinput_exp!=0) , ' ^{\circ}'*(preferredScale[whereParam[i_stats]]=='deg')) #Only writes out the power of ten if it is needed


                if write_tex_files:
                    outtext += '%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\\n'%(latexParamNames[whereParam[i_stats]],dist_type[i_path,i_stats], inside_sigma_text[i_path,i_stats] ,maxlike_text , mu_text , trueinput_text  , no_s_d_ , diff_sign_text[i_path,i_stats],sigma_percentage_text[i_path,i_stats])
             ############################
             ### Writing matrix files ###
             ############################
                
                #######
                ### matrix layout: mu, sigma, maximum likelihood, type of distribution (see note from 12/8 -15
                #######

                ### Saves file in the same directory as where the chains folder is found
                
            out_array = np.zeros([6,n_params])
            out_array[0] = stats_gaussian[i_path,:,0]
            out_array[1] = stats_gaussian[i_path,:,1]
            out_array[2] = stats_max_like[i_path]
            out_array[3] = type_of_distribution[i_path]
            out_array[4] = local_log_evidence[i_modes]
            out_array[5] = global_log_evidence
            np.savetxt('%s/%s_mode_%d.dat'%(stats_path_to_chains[i_path].split('chains/1-stats.dat')[0], output_filename , i_modes), out_array )
        if write_tex_files:
            outtext += '\\end{tabular}'
        
            write_output = open('%s.tex'%(output_filename),'w')
            write_output.write(outtext)
            write_output.close()



    #######################
    ### Analysing input ###
    #######################

    if read_output_matrix:
        print np.shape(input_data)
        print input_data_path
        print len(input_data_model_list)
        ### Degrees of Freedom for each type of coverage
        dict = {'full_log':162 , 'full_inconsistent':162 , 'optical':152 , 'xrays':12 , 'radio':138 , 'prebreak':152 , 'postbreak':152}
        
        true_input_list = [[51.,53.],[-2.,0.],[5.,15.],[1.,7.]]

        ########################
        ### Analysing biases ###
        ########################

        bias_after_model = np.zeros([7,32])

        for i_model in range(2):
            for j_model in range(2):
                for k_model in range(2):
                    for l_model in range(2):
                        for s_model in [0,2]:
                            ### Order in input_data: mu, sigma, maximum likelihood, type of distribution, local log-evidence (central value and 1 sigma), global log-evidence (central value and 1 sigma)
                            
                            blabla
                            ###################
####################3

#CONTINUE HERE!!!!
        
        


###############################################
#         Beginning of programme              #
###############################################  

#Importing user options and constants from options.py
inputData,fileExtention,breakName,inputTimeFormat,inputFreqUnit,plotInput,GRBlabel,GRBinput,allowDataRemoval,gridReducer,surfacePlot,thermalComp,photosphere,numericalOpticalDepth,analyticRhoProfile,useEATS,plotOutput,gridStart,gridEnd,gridStep,printProcess,daysOrSec,fluxAxis,reverseShock,opticalDepth,sanityCheckEATS,printStats,runOption,mockDistribution,numberOfMock,mockIntervalTime,mockIntervalFreq,createMock,gaussianError,gaussianOffset,offsetType,useData,plotComponents,plotMock,mockDim,fixedRSFSratio,figTypes,freqGrid,tgrid,cb,paramLimits,parametrar,allowPrint,printCrucial,analyseChains,livePoints,chi2_type = options.userOptions()

constants = np.array(options.inputConstants())    #If this statement is False, it means the user has choosen the parameters.txt file to be loaded

if allowPrint: print "Program started %dh %dm %ds"%(tidPre[3:6])

n_params = int(np.sum(parametrar))
currentPath = os.path.abspath('.')
abovePath = os.path.abspath('../')
homePath = os.path.expanduser('~')
whereParam = np.where(parametrar)[0]


surfRingsOut = 10
numRuns = 0
loadInputConstants = True #When this is true, the options.py file is loaded for constants. If not, the parameters.txt file is loaded
printPlot = False
plot_area = False
paramNames = ['epsilon','log10(epsilone_FS)','log10(epsilonp_FS)','log10(epsilone_RS)','log10(epsilonp_RS)','log10(E0)','log10(n)','s','Gamma0','log10(epsilonB_FS)','log10(epsilonB_RS)','p_FS','t0','theta0','alpha','tN','log(r0)','TMax','log(N0)','t_prompt','p_RS','R0 (cocoon)','N (cocoon)','t_outflow (cocoon)','theta0 (cocoon)','Gamma0 (cocoon)','z','WV']
latexParamNames = [r'$\epsilon$',r'$\log_{10}(\epsilon_{\rm e})$',r'$\log_{10}(\epsilon_{\rm p})$',r'$\log_{10}(\epsilon_{\rm e,RS})$',r'$\log_{10}(\epsilon_{\rm p,RS})$',r'$\log_{10}(E_0)$',r'$\log_{10}(n)$',r'$s$',r'$\log_{10}(\Gamma_0)$',r'$\log_{10}(\epsilon_{\rm B})$',r'$\log_{10}(\epsilon_{\rm B,RS})$',r'$p_{\rm FS}$',r'$t_0$',r'$\theta_0$',r'$\alpha$',r'$t_{\rm N}$',r'$\log_{10}(r_0)$',r'$T_{\rm Max}$',r'$\log_{10}(N_0)$',r'$\Delta t_{\rm prompt}$',r'$p_{\rm RS}$',r'$\log_{10}(R_{\rm 0,coc})$',r'$\log_{10}(N_{\rm coc})$',r'$\delta t_{\rm coc}$',r'$\theta_{\rm 0,coc}$',r'$\Gamma_{\rm 0,coc}$',r'$z$','WV']
latexParamNamesLin = [r'$\epsilon$',r'$\epsilon_{\rm e}$',r'$\epsilon_{\rm p}$',r'$\epsilon_{\rm e,RS}$',r'$\epsilon_{\rm p,RS}$',r'$E_{\rm 0,iso}$',r'$n$',r'$s$',r'$\Gamma_0$',r'$\epsilon_{\rm B}$',r'$\epsilon_{\rm B,RS}$',r'$p_{\rm FS}$',r'$t_0$',r'$\theta_0 (^{\circ})$',r'$\alpha (^{\circ})$',r'$t_{\rm N}$',r'$r_0$',r'$T_{\rm Max}$',r'$N_0$',r'$\Delta t_{\rm prompt}$',r'$p_{\rm RS}$',r'$R_{\rm 0,coc}$',r'$N_{\rm coc}$',r'$\delta t_{\rm coc}$',r'$\theta_{\rm 0,coc} (^{\circ})$',r'$\Gamma_{\rm 0,coc}$',r'$z$','WV']
paramNamesShort = ['epsilon','epsilon_e','epsilon_p','epsilon_e_RS','epsilon_p_RS','E0','n','s','Gamma0','epsilon_B','epsilon_B_RS','p','t0','theta0','alpha','tN','r0','Tmax','N0','t_prompt','p_RS','R0_coc','N_coc','t_outflow','theta0_coc','Gamma0_coc','z','WV']
preferredScale = ['log','log','log','log','log','log','log','lin','log','log','log','lin','None','deg','deg','None','None','None','None','log','lin','log','log','log','deg','log','lin','lin']
constantsLength = len(paramNames) #The number of elements in contants list



### Handling input arguments ###


nArgs = len(sys.argv)
FS_only , RS_only , thermal_only = False,False,False
bin_data_bool = False
if nArgs > 1: 
    inArg = sys.argv
    
    #Help input argument
    if (inArg[1] == "h") | (inArg[1] == "help"):
        print "This is the help section, called by commandline arguments h or help.\n\n--------------------------------\nAvailable commandline arguments: \n\nburst= Type the date of the GRB\noption= \n    LC: lightcurve production\n    fit: run the fitting routine\n    marginal: plot 2D marginal plots\n    plot-marginal: advanced marginal plotter. Lets the user plot many posteriors on top of each other\n    print-stats: Produce a LaTeX table with stats from MultiNest output\n    read-stats: Produce a LaTeX table with the properties of the probability distributions. LaTeX output requires package \'color\'\n    the default option is set in options.py, named runOption=''\n\nload=\n    output: load constants from the parameters.txt file, printed by the fitting routine\n    chains: Load best-fit parameters from MultiNest output\n    Default is to load the parameters from the options.py file\n\nexcept=\n    FS: Plot the forward-shock only\n    RS: Plot the reverse-shock component only\n    thermal: Plot the thermal component only\n    These above may be combined by calling except= many times. Option expept is best suited when printing lightcurve output files (print-write). Save each component in a specific file, and load all files when running print-read. See manual for mor details, by entering commandline argument man. This option is still under construction\n\nbin=\n    Enter number of data points per bin. Default is no binning.\n\nband=\n    xrays/x-rays: Load X-ray lightcurves >0.1 keV (2.4e16 Hz)\n    UV: Load UV lightcurves 3.2 eV - 0.1 keV (7.7e14 - 2.4e16 Hz)\n    UV-optical/optical-UV: Load optical and UV lightcurves 1.6eV - 3.2eV\n    optical: Load optical lightcurves 0.4eV - 4eV\n    radio-mm: Load radio-mm lightcurves < 1.6eV\n    radio: Loading radio lightcurves (<600 GHz)\n\nProgramme is called by typing \'(i)python mainFit.py ->commandline arguments<-\'\n\nFor other questions, mail Andreas Johansson on johansson.mcquack@gmail.com\nNow exiting.  "
        raise SystemExit(0)
    elif inArg[1] == 'man':
        print 'There is no manual pages writted yet. Please contact Andreas Johansson via johansson.mcquack@gmail.com\nNow exiting.'
        raise SystemExit(0)
    #Other input arguments
    for i in range(1,nArgs):
        inputArg = inArg[i].split('=')
        #Setting burst name

        if inputArg[0] == 'except':
            if inputArg[1] == 'FS': FS_only = True
            elif inputArg[1] == 'RS': RS_only = True
            elif inputArg[1] == 'thermal': thermal_only = True
            else: print 'No except I can understand!\n\nPossible expects:\n\'FS\' - plot forward shock only\n\'RS\' - plot reverse shock only'

        elif inputArg[0] == 'bin':
            try: 
                bin_data = int(inputArg[1])
                bin_data_bool = True  #Boolean variable to determine whether to bin input data or not.
            except: 
                print 'Invalid number of points per bin. Type bin=number_of_points_per_bin. Now exiting'
                raise SystemExit(0)
        elif inputArg[0] == 'burst':
            GRBlabel = inputArg[1]
        #Setting runOption
        elif inputArg[0] == 'option':
            if inputArg[1] == 'print-write':
                runOption = 'LC'
                printPlot = True
                plotComponents = False
            elif inputArg[1] == 'print-read':
                runOption = 'LC'
                plotComponents = True
                printPlot = True
            elif inputArg[1] == 'area':
                runOption = 'LC'
                plot_area = True
            else:
                runOption = inputArg[1]
 
        elif inputArg[0] == 'load':
            print "Loading parameters.txt"
            if inputArg[1] == 'output': #Loading the output file parameters.txt created in the lowest-chi2 routine
                try:
                    os.path.isfile('parameters.txt')
                    ladda = open('parameters.txt')
                    loadInText = ladda.read().split('\n')
                    ladda.close()
#                    constants = np.zeros(constantsLength)
                    for j in range(constantsLength):
                        constants[j] = float(loadInText[j].split('=')[1])
                    surfRingsOut = int(loadInText[constantsLength].split('=')[1])
                    loadInputConstants = False
                except: print "Could not find file parameters.txt. Loading the options.py file instead"
            elif inputArg[1] == 'chains':
                ### Reading best-fit values from chains/1-stats.dat
                print 'Loading best-fit parameters from MultiNest output'
                try:
                    open_stats = open('chains/1-stats.dat')
                    read_stats = open_stats.read().split('Dim No.')[2].split('\n')
                    open_stats.close()
                    for i_chains in range(1,n_params+1):
                        constants[whereParam[i_chains-1]] = float(read_stats[i_chains].split()[1])
                except:
                    print 'Failed when reading chains/1-stats.dat. Either file does not exist, or the layout of the MultiNest output has been changed. Now exiting.'
                    raise SystemExit(0)
            else: print "WARNING! No argument %s"%inputArg[1]
        elif inputArg[0] == 'band':
            bandLimit = np.zeros(2)   #bandLimit[0] lower limit, bandLimit[1] upper limit
            if (inputArg[1] == 'x-ray') or (inputArg[1] == 'X-ray') or (inputArg[1] == 'x-rays') or (inputArg[1] == 'X-rays') or (inputArg[1] == 'xrays'):
                print "Loading X-ray lightcurves >0.1 keV (2.4e16 Hz)"
                bandLimit[0] = 2.4e16
                bandLimit[1] = 0
                bandLabel = "X-rays >0.1 keV (>2.4e7 GHz)"
            elif (inputArg[1] == 'UV') or (inputArg[1] == 'uv'):
                print "Loading UV lightcurves 3.2 eV - 0.1 keV (7.7e14 - 2.4e16 Hz)"
                bandLimit[0] = 7.7e14
                bandLimit[1] = 2.4e16
                bandLabel = "UV 3.2-100 eV (7.7e5 - 2.4e7 GHz)"
            elif (inputArg[1] == 'optical-UV') or (inputArg[1] == 'UV-optical'):
                print "Loading optical and UV lightcurves 1.6eV - 3.2eV"
                bandLimit[0] = 2.e14
                bandLimit[1] = 2.4e16
                bandLabel = "UV and optical 1.6-100 eV (2e5 - 2.4e7 GHz)"
            elif (inputArg[1] == 'optical'):
                print "Loading optical lightcurves 0.4eV - 4eV"
                bandLimit[0] = 1.e14
                bandLimit[1] = 1.e15
                bandLabel = "optical 0.4-4 eV (1e5/3cm - 1e6/300nm GHz)"
            elif (inputArg[1] == 'radio-mm'):
                print "Loading radio-mm lightcurves < 1.6eV"
                bandLimit[0] = 0
                bandLimit[1] = 2.e14
                bandLabel = "radio - submm - IR < 1.6 eV (<2.4e5 GHz)"
            elif (inputArg[1] == 'radio'):
                print "Loading radio lightcurves"
                bandLimit[0] = 0
                bandLimit[1] = 6e11 #lambda = 2mm
                bandlabel = 'radio - submm (<2mm ; <600 GH)'
            elif (len(inputArg[1].split('-'))==2):
                bandLimit[0] = float(inputArg[1].split('-')[0])
                bandLimit[1] = float(inputArg[1].split('-')[1])
                printBand1, printBand1e, printBand2 , printBand2e = bandLimit[0]/1e9 , 0 , bandLimit[1]/1e9 , 0
                while printBand1 >= 10:
                    printBand1 /= 10
                    printBand1e += 1
                while printBand2 >= 10:
                    printBand2 /= 10
                    printBand2e += 1
                print "Loading %se%d - %se%d GHz"%(printBand1,printBand1e,printBand2,printBand2e)
                bandLabel = "%se%d - %se%d GHz"%(printBand1,printBand1e,printBand2,printBand2e)
            else:
                try:
                    
                    bandLimit[0] = float(inputArg[1])
                    bandLimit[1] = float(inputArg[1])
                    print "Loading %s"%bandLimit[0]
                except:
                    print "Bad band input. Loading all lightcurves..."

    print_read = printPlot and plotComponents
    print_write = printPlot and (not plotComponents)


else:
    if allowPrint: print "No input arguments"


if (createMock) & (runOption == 'LC'): 
    useData = False    #Overrun
    plotMock = True
if (createMock == False) & (useData == False): plotMock = False   #Overrun
if runOption != 'fit': plotOutput == False
if (createMock) & (runOption != 'LC'): createMock = False



#Creating pathways
inputData = '%s%s%sfit/'%(inputData , '/'*(inputData.split('/')[-1]!='') , GRBlabel)
if not os.path.isdir(inputData):
    #if (runOption == 'LC') & createMock: 
    #    os.system('mkdir %s'%inputData)
    if runOption == 'fit': raw_input('No input data found. Input data should be saved in %s'%inputData)
elif createMock and (runOption=='LC'):
    deleteDir = raw_input('Save files in existing directory %s? ([y]/n)'%inputData)
    if deleteDir == '': deleteDir = 'y'
    if (deleteDir != 'y') and (deleteDir != 'Y'): 
        print "Relabel the burst in file option.py or by adding flag burst=\nNow exiting"
        raise SystemExit(0)
        

color_68='#616161'
color_95='#cfcfcf'


#mockPath = '%sMock/'%inputData
#if (runOption == 'LC') & createMock:
#    if os.path.isdir(inputData) == False: os.system('mkdir %s'%inputData)


if os.path.isfile('interactive_log.txt') & (runOption=='fit'): os.system('rm interactive_log.txt')
while allowDataRemoval & allowPrint:
    removePrevious = raw_input('OK removing previous data? [default=no, y=yes]: ') 
    if removePrevious == 'y':
        os.system('rm -r chains/') 
        useArchive = False #Whether files should be saved to archive after run or not. False: write to archive
        break
    elif removePrevious == '':
        useArchive = True
        break
    else: print "Bad input"

if os.path.isdir('chains/') == False: os.system('mkdir chains')
os.system('chmod a+rw chains/')


R = np.logspace(gridStart,gridEnd,gridStep)


parameterIndeces,parIndCount = np.array([0]*np.sum(parametrar)),0
colourCycle = ['b','g','r','c','m','y','k']     #Cycle of colours in plotting. Matplotlib standard cycle
scalePlotTime = {'d': 86400. , 'h' : 3600. , 'm' : 60. , 's' : 1. }
scaleFluxAxis = {'mJy' : 1.e3 , 'Jy' : 1. }

for parI in range(len(parametrar)):
    if parametrar[parI] == True: 
        parameterIndeces[parIndCount] = parI
        parIndCount += 1

if runOption == 'plot':
    plotter()
    raise SystemExit(0)
elif runOption == 'read-stats':
    gather_multinest_stats()
    raise SystemExit(0)
elif runOption == 'plot-marginal':
    from plotMarginal import plot_contours
    from plotMarginal import plotMarginal
    from matplotlib import pyplot as plt
#    from plot_marg_txt import plot_marg_txt
    choise = input('Plot marginal plots in a gathered sheet (1) or separately with contours (2)?: ')
    if choise == 1:
        plotMarginal(n_params)
    elif choise == 2: 
#        plot_saved = raw_input('Plot saved txt files? ([y]/n): ')
#        plot_saved = (plot_saved == '') or (plot_saved == 'Y') or (plot_saved == 'y')
#        if plot_saved:
#            plot_marg_txt(paramNamesShort,whereParam)
#            raise SystemExit(0)
        
        number_of_plots = input('Number of plots: ')
        path_to_chains = np.zeros(number_of_plots,dtype='S256')
        line_styles = np.zeros(number_of_plots,dtype='S5')
        for set_path in range(number_of_plots):
            path_to_chains[set_path] = raw_input('Path to %d%s chains folder: '%(set_path+1,'st'*(set_path==0)+'nd'*(set_path==1)+'rd'*(set_path==2)+'th'*(set_path>2)))
            line_styles[set_path] = raw_input('Line style [-]: ')
            if line_styles[set_path] == '': line_styles[set_path] = '-'

        if not os.path.isdir('Figures'): 
            os.system('mkdir Figures')
            print 'Created directory Figures'
        marginal_prefix = raw_input('Choose name prefix for these files: ')
        if not os.path.isdir('Figures/%s'%marginal_prefix): 
            os.system('mkdir Figures/%s'%marginal_prefix)
            print 'Created directory Figures/%s'%marginal_prefix

   ### Binning and plotting. Change manually in code to alternate the type of plot output                                                                                       

        n_bins = 20#int(np.sqrt(n_points)) #Number of bins on the axis of each parameter     
        cons_plot_num = 1 #Number to assign different plot number to every plot
        plot_num = np.zeros([n_params,n_params],dtype=int)

        for plot_path in range(number_of_plots):
            probability , parameters = plot_contours(n_params,os.path.abspath(path_to_chains[plot_path]))

            n_modes = np.shape(probability)[0]
            density = np.zeros([n_modes,n_params,n_params,n_bins,n_bins])
        
            for i_modes in range(n_modes):


                for rescale_par in range(n_params):
                    if preferredScale[whereParam[rescale_par]] == 'deg':
                        parameters[i_modes,rescale_par] *= 180 / np.pi
                    #elif preferredScale[whereParam[rescale_par]] == 'log':
                    #    parameters[i_modes,rescale_par] = 10**parameters[i_modes,rescale_par]
                sum_prob = 0.
                prob_temp = np.copy(probability[i_modes])


                ### The arrays parameters and probability may contain zeros (if there are multiple modes), and low-probability outliers. We pick out the points within two sigma, and use them as the bin limits
                while sum_prob < 0.95:
                    prob_max = np.argmax(prob_temp)
                    sum_prob += prob_temp[prob_max]
                    prob_temp[prob_max] -= 1
                    prob_ind = np.where(prob_temp < 0)
            

                for x in range(n_params):
                    for y in range(n_params):
                        if x >= y: continue #We do not want to plot the parameter on itself                                                                                              
                        if plot_path == 0:
                            plot_num[x,y] = cons_plot_num
                        exec('fig_%d_%d = plt.figure(%d)'%(x,y,plot_num[x,y]))
#                        if preferredScale[whereParam[x]] == 'log':
#                            bin_limits_x = np.logspace(min(np.log10(parameters[i_modes,x][prob_ind])),max(np.log10(parameters[i_modes,x][prob_ind])),n_bins+1)
#                            print bin_limits
#                            raw_input(paramNames[whereParam[x]])
#                        else:
                        bin_limits_x = np.linspace(min(parameters[i_modes,x][prob_ind]),max(parameters[i_modes,x][prob_ind]),n_bins+1)
#                        if preferredScale[whereParam[y]] == 'log':
#                            bin_limits_y = np.logspace(min(np.log10(parameters[i_modes,y][prob_ind])),max(np.log10(parameters[i_modes,y][prob_ind])),n_bins+1)
#                        else:
                        bin_limits_y = np.linspace(min(parameters[i_modes,y][prob_ind]),max(parameters[i_modes,y][prob_ind]),n_bins+1)
                        
                        
                        ### Binning

                        for x_bins in range(n_bins):
                            for y_bins in range(n_bins):
                                density[i_modes,x,y,x_bins,y_bins] = np.sum(probability[i_modes,np.where((parameters[i_modes,x]>bin_limits_x[x_bins]) & (parameters[i_modes,x]<bin_limits_x[x_bins+1]) & ((parameters[i_modes,y]>bin_limits_y[y_bins])) & ((parameters[i_modes,y]<bin_limits_y[y_bins+1])))])
                        plot_x , plot_y = np.meshgrid(np.linspace(bin_limits_x[0],bin_limits_x[-1],n_bins) , np.linspace(bin_limits_y[0],bin_limits_y[-1],n_bins))

#                        plt.subplot(n_params,n_params,y+x*(n_params+1)+1)
                        exec('fig_%d_%d'%(x,y))
                        
                        if preferredScale[whereParam[x]] == 'log':
                            plot_x = 10**plot_x
                            plt.xscale('log')
                        if preferredScale[whereParam[y]] == 'log':
                            plot_y = 10**plot_y
                            plt.yscale('log')
                        if paramNames[whereParam[x]] == 'Gamma0 (cocoon)':
                            plot_x = np.sqrt(plot_x**2 - 1)
                        elif paramNames[whereParam[y]] == 'Gamma0 (cocoon)':
                            plot_y = np.sqrt(plot_y**2 - 1)
                        plt.contour(plot_x,plot_y, density[i_modes,x,y],10,linestyles=line_styles[plot_path])
                        plt.xlabel(latexParamNamesLin[whereParam[x]])
                        plt.ylabel(latexParamNamesLin[whereParam[y]])
                        cons_plot_num += 1
                        #color_cplot = plt.colorbar(cplot,shrink=0.8,extend='both')
                        
                        np.savetxt('Figures/%s_%s_xaxis_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , plot_x)
                        np.savetxt('Figures/%s_%s_yaxis_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , plot_y)
                        np.savetxt('Figures/%s_%s_density_%d.txt'%(paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]] , plot_path) , density[i_modes,x,y])
                        
                        if plot_path == (number_of_plots-1):
                            exec('fig_%d_%d.savefig(\'Figures/%s/%s_%s.%s\')'%(x,y,marginal_prefix,paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]],figTypes))
                            print 'Saved figure Figures/%s/%s_%s.%s'%(marginal_prefix,paramNamesShort[whereParam[x]],paramNamesShort[whereParam[y]],figTypes)
        #plt.show()
        
        print np.sum(density)



    else:
        print 'Invalid option %d. Now exiting.'%choise
    raise SystemExit(0)


### Load data to plot uncertainty range
def load_sigma(sigma_prob_range,sigma_prob_input,load_area_name):
    for sigma_prob in sigma_prob_range:
        load_name_lowest = '%s/lowest_%d.txt'%(load_area_name,int(100*sigma_prob))
        load_name_highest = '%s/highest_%d.txt'%(load_area_name,int(100*sigma_prob))
        load_name_tempgrid = '%s/temp_grid_%d.txt'%(load_area_name,int(100*sigma_prob))
        if sigma_prob == 0.68:
            lowest_LC_68 = np.loadtxt(load_name_lowest)
            highest_LC_68 = np.loadtxt(load_name_highest)
            sigmaTempGrid_68 = np.loadtxt(load_name_tempgrid)
        elif sigma_prob == 0.95: 
            lowest_LC_95 = np.loadtxt(load_name_lowest)
            highest_LC_95 = np.loadtxt(load_name_highest)
            sigmaTempGrid_95 = np.loadtxt(load_name_tempgrid)
    if sigma_prob_input == 3: return lowest_LC_68, highest_LC_68, sigmaTempGrid_68, lowest_LC_95, highest_LC_95, sigmaTempGrid_95 #Plot both 68% and 95 %
    elif sigma_prob_input == 1: return lowest_LC_68, highest_LC_68, sigmaTempGrid_68,[],[],[]  #Plot only 68%
    else: return [],[],[], lowest_LC_95, highest_LC_95, sigmaTempGrid_95 #Plot only 95%


#Loading data. If we are creating mock observations, data will not be loaded
if ((runOption != 'LC') | (createMock == False)) and (runOption != 'plot'):
    if (createMock == False) & ((runOption == 'LC') & (useData == False)):  #If we only want to plot a certain model
        if mockDistribution == 'log': tdata = 10**np.array([np.linspace(np.log10(mockIntervalTime[0]),np.log10(mockIntervalTime[1]),numberOfMock)]*len(freqGrid)) 
        elif mockDistribution == 'lin': tdata = np.array([np.linspace(mockIntervalTime[0],mockIntervalTime[1],numberOfMock)]*len(freqGrid)) 
        else: 
            try:tdata = np.array[mockDistribution]
            except: 
                if allowPrint: print "Bad mockDistribution input. Exiting."
                raise SystemExit(0)
        freqInput =  np.array(freqGrid)
    else:tdata,FdataInput,errorbarInput,errorIsArr,freqInput = loadData()

    if createMock: freq = freqMock
        
            
    else: freq = freqInput
    print "Number of data points = %d"%numberOfPoints
        

#When wished output is in the temporal grid, the frequency might be a scalar, and the same for the time when the frequency is resolved. To keep the code general, we now make sure that these quantities alwas are arrays
if mockDim == 'E': 
    if createMock: tdata = np.array(tgrid)
    try: iterationLength = len(tdata)    #Makes sure that input tdata is an array
    except: 
        tdata = np.array([tdata])
        iterationLength = 1
    numberOfEmpties = np.zeros(iterationLength)
    if createMock == False:
        for i in range(iterationLength):
            numberOfEmpties[i] = np.count_nonzero(freq[i]+1)
elif mockDim == 'T':
    if createMock: freq = np.array(freqGrid)
    try: iterationLength = len(freq)
    except: 
        freq = np.array([freq])  #Same as above
        iterationLength = 1
    numberOfEmpties = np.zeros(iterationLength, dtype=int)
    if createMock == False:
        for i in range(iterationLength):
            numberOfEmpties[i] = np.count_nonzero(tdata[i]+1.)



### Plotting the 1-sigma area around the lightcurves

def plot_area_func():
    
    load_area_name = raw_input('Type GRB label to store uncertainty area files. If files exist for this name, those will be loaded. Leave blank to use the input burst name: ') #This input is used when plotting the areas
    ### Picking what uncertainty range to plot
    if load_area_name == '': load_area_name = GRBlabel
    try:sigma_prob_input = input('Plot 1-sigma ([1]) or 1-sigma and 2-sigma(2)?: ')
    except: #Input must be an integer
        sigma_prob_input = 1
    if sigma_prob_input == 1: sigma_prob_range = [0.68]
    elif sigma_prob_input == 2: sigma_prob_range = [0.68,0.95]
    else: #If bad input
        print 'Bad input %s! Now exiting'%sigma_prob_input
        raise SystemExit(0)
    if os.path.isdir(load_area_name):
        lowest_LC_68, highest_LC_68, sigmaTempGrid_68, lowest_LC_95, highest_LC_95, sigmaTempGrid_95 = load_sigma(sigma_prob_range,sigma_prob_input,load_area_name)
    else:
        for sigma_prob in sigma_prob_range:
            all_sigma_points = get_one_sigma(n_params,0.68)
        
        ### Picking 10% of the points at random

            sigma_points = all_sigma_points[random.sample(range(len(all_sigma_points)),len(all_sigma_points)/10)]
            print 'Getting %d points from data set, using %d points for plotting uncertainty range'%(len(all_sigma_points),len(sigma_points))

            for sigma_LC in range(len(sigma_points)):

        ### Setting constants from the one-sigma points
                
                constants[whereParam] = sigma_points[sigma_LC,2:]
                    
                lightcurve_sigma,hej,hej,sigmaTempGrid = modelFunc(R,constants,fixedRSFSratio,reverseShock,'one-sigma',useEATS,thermalComp,photosphere,printStats,allowPrint,tdata,FdataInput,errorbarInput,freq,iterationLength,numberOfEmpties,createMock,opticalDepth,plotComponents,daysOrSec,fluxAxis,numberOfPoints,printProcess,surfRingsOut,np.sum(parametrar),FS_only,RS_only,thermal_only,chi2_type)
                    
                    
                if sigma_LC == 0:
                    lowest_LC = np.copy(lightcurve_sigma)
                    highest_LC = np.copy(lightcurve_sigma)
                else:                

                ### Assigning values if extremes are reached
                    for extremes in range(len(lightcurve_sigma)):
                        where_lower = np.where(lightcurve_sigma[extremes] < lowest_LC[extremes])[0]
                        where_higher = np.where(lightcurve_sigma[extremes] > highest_LC[extremes])[0]
                        lowest_LC[extremes,where_lower] = lightcurve_sigma[extremes,where_lower]
                        highest_LC[extremes,where_higher] = lightcurve_sigma[extremes,where_higher]
                    print sigma_LC
            
            if not os.path.isdir(load_area_name): os.system('mkdir %s/'%load_area_name)
            lowest_filename = '%s/lowest_%d.txt'%(load_area_name,int(100*sigma_prob))
            highest_filename = '%s/highest_%d.txt'%(load_area_name,int(100*sigma_prob))
            tempgrid_filename = '%s/temp_grid_%d.txt'%(load_area_name,int(100*sigma_prob))
            np.savetxt(lowest_filename,lowest_LC)
            np.savetxt(highest_filename,highest_LC)
            np.savetxt(tempgrid_filename,sigmaTempGrid)
            print 'Files %s, %s and %s are written'%(lowest_filename,highest_filename,tempgrid_filename)
                
            if sigma_prob == 0.68:
                lowest_LC_68 = np.copy(lowest_LC)
                highest_LC_68 = np.copy(highest_LC)
                sigmaTempGrid_68 = np.copy(sigmaTempGrid)
            elif sigma_prob == 0.95:
                lowest_LC_95 = np.copy(lowest_LC)
                highest_LC_95 = np.copy(highest_LC)
                sigmaTempGrid_95 = np.copy(sigmaTempGrid)
                    
        np.savetxt('%s/freq.txt'%load_area_name,freq)

    return sigma_prob_range,sigma_prob_input,load_area_name



def lightcurve_production(freq,numberOfMock,numberOfEmpties,FdataInput,tdata,errorbarInput):
    from matplotlib import pyplot as plt
    from matplotlib import gridspec
    import matplotlib as mpl
    from plot_sigma import get_one_sigma
    from plotMarginal import plotMarginal
#    global FdataInput,errorbarInput,numberOfPoints,tdata
    ccgs = 2.99792458e10   #Speed of light in CGS units
    if loadInputConstants: surfRingsOut = 200   #To assure we have enough with rings in EATS integrator
    if createMock:
        FdataInput,errorbarInput,numberOfPoints = [],[],0
    #Creating temporal grid...
        if mockDistribution == 'log': #Logarithmic distance between mock data points
            if mockDim == 'T':     tdata = np.array([10**np.linspace(np.log10(mockIntervalTime[0]),np.log10(mockIntervalTime[1]),numberOfMock)]*len(freq))
            elif mockDim == 'E':   freq = np.array([10**np.linspace(np.log10(mockIntervalFreq[0]),np.log10(mockIntervalFreq[1]),numberOfMock)]*len(tdata))
        elif mockDistribution == 'lin':
            if mockDim == 'T':     tdata = np.array([np.linspace(mockIntervalTime[0],mockIntervalTime[1],numberOfMock)]*len(freq))
            elif mockDim == 'E':   freq = np.array([np.linspace(mockIntervalFreq[0],mockIntervalFreq[1],numberOfMock)]*len(tdata))
        else: 
            if mockDim == 'T':     tdata = mockDistribution
            elif mockDim == 'E':   freq = mockDistribution
            if nl > 1: numberOfMock = len(mockDistribution[0])
            else: numberOfMock = len(mockDistribution)
        numberOfEmpties = np.array([numberOfMock]*iterationLength)

        
    nl = len(freq)

    
    


    ### Plotting the lightcurves

    if not print_read:    # If this is false, the user want's to read the printed data and plot it
        if createMock:lightcurve , tempGrid = logLikelihood([],[],[],FdataInput,tdata,errorbarInput,numberOfEmpties)
        else: lightcurve , tempGrid = logLikelihood([],[],[])

    
    ### Create and save mock observations

    if createMock:
        if mockDim == 'T':    mockGrid,mockLC = tdata,lightcurve   #In janskys
        elif mockDim == 'E':  mockGrid,mockLC = freq,lightcurve
        FdataInput = lightcurve
        mockError = mockLC * gaussianError      #For now
        #Gaussian offset
        if gaussianOffset: 
            store_seed = np.zeros([len(mockLC),numberOfMock])
            for i in range(len(mockLC)):
                
                for j in range(numberOfMock):
                    #Set seed
                    time_now = time.time()
                    microseconds = int((time_now - int(time_now)) * 1e6)
                    store_seed[i,j] = microseconds
                    random.seed(microseconds)
                    #Create data points
                    if offsetType == 'lin':
                        mockLC[i,j] = random.gauss(mockLC[i,j],mockLC[i,j]*gaussianError)
                        while mockLC[i,j] <= 0: mockLC[i,j] = random.gauss(mockLC[i,j],mockLC[i,j]*gaussianError) #Preventing negative flux
                    elif offsetType == 'log':                        
                        mockLC[i,j] = np.random.lognormal(np.log(mockLC[i,j]),np.e**mockLC[i,j]*gaussianError)
                        while mockLC[i,j] <= 0: mockLC[i,j] = random.gauss(np.log10(mockLC[i,j]),np.log10(np.abs(mockLC[i,j]))*gaussianError) #Preventing negative flux
                    else:
                        print 'Bad offset type \'%s\'. Now exiting'
                        raise SystemError(0)
            
        print os.path.isdir(inputData)
        print inputData
        if os.path.isdir(inputData):
            loopAgain = True
            while loopAgain:
                optIn = raw_input('Erase files in %s before saving files there? ([y]/n): '%inputData)
                if optIn == '': optIn = 'y'
                if (optIn == 'y') | (optIn == 'Y'):
                    os.system('rm -r %s*'%inputData)
                    loopAgain = False
                elif (optIn == 'n') | (optIn == 'N'): loopAgain = False
        for iNu in range(iterationLength):
            outText = "%d\n"%(numberOfMock-(len(mockLC[iNu]) - np.count_nonzero(mockLC[iNu])))
            for oi in range(np.count_nonzero(mockLC[iNu]+1)): 
                #if mockLC[iNu,oi] != 0.:outText += "%s %s %s\n"%(mockGrid[iNu,oi],mockLC[iNu,oi],mockLC[iNu,oi]*0.1)#mockError[iNu,oi]) #Producing file with a fix errorbar (not necessarily the correct errorbar)
                if mockLC[iNu,oi] != 0.:outText += "%s %s %s\n"%(mockGrid[iNu,oi],mockLC[iNu,oi],mockError[iNu,oi])
            if mockDim == 'T':
                if (os.path.isdir('%s'%inputData)==False): os.system('mkdir %s'%inputData)
                if os.path.isdir('%s'%inputData) == False: raw_input('Couldn\'nt create directory %s%stimeResolved! Create it manually before continuing!'%(dataFiles,'/'*(dataFiles.split('/')[-1] != '')  ) )
                os.system('touch %s%slc.%s'%(inputData,freq[iNu],fileExtention))
                a = open('%s%slc.%s'%(inputData,freq[iNu],fileExtention),'w')
            elif mockDim == 'E': 
                if os.path.isdir('%s%sfrequencyResolved/'%(inputData, '/'*(inputData.split('/')[-1] != '') )) == False: os.system('mkdir %s%sfrequencyResolved'%(inputData,'/'*(inputData.split('/')[-1] != '') ))
                a = open('%s%sfrequencyResolved/%slc.%s'%(inputData,'/'*(inputData.split('/')[-1] != ''),tdata[iNu],fileExtention),'w')
            a.write(outText)
            a.close()
        np.savetxt('%s/seed.txt'%inputData,store_seed)  #Storing seed used to produce mock data
        print "Mock data produced. Now exiting"
        raise SystemExit(0)
       

    #####################
    ### Plotting data ###
    #####################

    if numberOfMock == 1: plotDims = raw_input('Do you want to plot the lightcurve with frequency on the x-axis? (y/n): ')
    else: plotDims = 'n'
            
     ### Writing print-outs. If user input option=print-write  -  print plot to a txt file

    if print_read:  #Loading print-outs
        fileNameIn = raw_input('Type file name or list of file names, separated with a space [print_plot.txt]: ')
        if fileNameIn == '': fileNameIn = 'print_plot.txt'
        fileName = fileNameIn.split()
    elif print_write: #Writing print-outs
        fileNameIn = raw_input('Enter name of file to save plot data in [default=print_plot.txt]: ')
        if fileNameIn == '': fileName = 'print_plot.txt'
        else: fileName = fileNameIn

        if os.path.isfile(fileName):
            print "File %s exists. Adding plot data after existing data."%fileName
            readPlotPrint = open(fileName)
            readPrint = readPlotPrint.read()
            readPlotPrint.close()
        else: 
            readPrint = ''
            os.system('touch %s'%fileName)


    ### Plotting print-outs ###

    if not print_write: #If the option is not to write a print-out file
        if print_read: #Read the print-out
            
            bandNames = ['U-band','B-band','V-band','G-band','R-band','I-band','Z-band','Y-band','J-band','H-band','K-band','L-band','M-band','N-band','Q-band']
            bandWaveLengths = [365.,445.,551.,605.,658.,806.,900.,1020.,1220.,1630.,2190.,3450.,4750.,10500.,2100.]
            use_standard_colors_input = raw_input('Use standard color scheme? ([y]/n): ')
            if (use_standard_colors_input=='y') or (use_standard_colors_input=='Y') or (use_standard_colors_input==''): use_standard_colors = True
            else: use_standard_colors = False
            bandStandardColors = ['k','c','m','g','b','y','r','c','k','y','r','g','m','b','c','k']
            standardFreqs = [1.9e9,4.8e9,8.4e9,22.5e9,43e9,100e9,300e9,1.4e14,2.5e14,3.7e14,4.6e14,4.8e14,5.4e14,6.7e14,8.2e14,2.4e18]
            legend = ['']*iterationLength

            lineFreq = np.array([])  # Each frequency in print-out file gets an entry in this array
            lineColorMem = np.array([])
            scaleLC = np.array([],dtype=float)
            lineTypeMem = np.array([],dtype=float)
            #legend = np.array([])
            bandName = np.array([],dtype=str)
            freqMem = np.array([[]])
            tempMem = np.array([[],[]])
            fluxMem = np.array([])
            tempMemIndex = np.array([0],dtype=int)

            
            plot_area_input = raw_input('Plot uncertainty areas? ([y]/n): ')
            if (plot_area_input == 'y') or (plot_area_input == 'Y') or (plot_area_input == ''):
                sigma_prob_range,sigma_prob_input,load_area_name = plot_area_func()
                
                lowest_LC_68, highest_LC_68, sigmaTempGrid_68, lowest_LC_95, highest_LC_95, sigmaTempGrid_95 = load_sigma(sigma_prob_range,sigma_prob_input,load_area_name)
                if os.path.isfile('%s/freq.txt'%load_area_name): sigma_freq = np.loadtxt('%s/freq.txt'%load_area_name)
                else: sigma_freq = np.loadtxt('%s/freq.txt'%load_area_name)

                plot_area = True
            else:
                plot_area = False

            scaleLightcurves = raw_input("Scale the lightcurves? ([y]/n): ")
            if scaleLightcurves == '': scaleLightcurves = 'y'
            if (scaleLightcurves == 'y') or (scaleLightcurves == 'Y'):
                ### Loading data file with scale factors for lightcurve scaling
                loadScaleData = raw_input('Name of list file with scale input. Leave blank for manual input: ')
                if loadScaleData != '':
                    try:
                        openScaleData = open(loadScaleData)
                        scaleDataInput = openScaleData.read().split('\n')
                        openScaleData.close()
                    
                        ### Reading in scale data to variable scaleData (nx2 matrix)
                        scaleData = np.zeros([len(scaleDataInput)-1,2],dtype=float)
                        for iScale in range(len(scaleData)):
                            scaleData[iScale] = map(float,scaleDataInput[iScale].split())
                    except:
                        print 'No such file. Now exiting.'
                        raise SystemExit(0)
                else: #If input is left empty, code will assume manual input and giving the opportunity to save input to a text file for later use
                    saveScaleData = ''
            for iFile in range(len(fileName)):
                readPrintOpen = open(fileName[iFile])
                readPrint = readPrintOpen.read().split('\n')
                readPrintOpen.close()
                        
                readLength = len(readPrint) / 3
                readLineLength = len(readPrint[1].split())
                tempGrid,lightcurve,readFreq,linecolor = np.zeros([readLength,readLineLength]),np.zeros([readLength,readLineLength]),np.zeros(readLength),np.zeros(readLength,dtype=str)
                
                lineType = raw_input('Choose line type for file %s (-,--,:,-.): '%fileName[iFile])
                while not ((lineType == '-') | (lineType == '-.') | (lineType == ':') | (lineType == '--')):
                    lineType = raw_input("Bad line type \'%s\'! Try again: "%lineType)

                for iLine in range(0,len(readPrint),3):  #Looping over input frequencies
                    if readPrint[iLine] == '': break
                    readFreq[iLine/3] = float(readPrint[iLine].split(':')[1])
                    readoutFreq , readFreqExp = readFreq[iLine/3] , 0

                    while readoutFreq > 10: #Reducing frequency to a reader-friendly form
                        readoutFreq /= 10
                        readFreqExp += 1

                    tempGrid[iLine/3] = map(float,readPrint[iLine+1].split())
                    lightcurve[iLine/3] = map(float,readPrint[iLine+2].split())
                    lineTypeMem = np.append(lineTypeMem,lineType)

                    
                    tempMem = np.append(tempMem,tempGrid[iLine/3])
                    fluxMem = np.append([fluxMem],[[lightcurve[iLine/3]]])
                    freqMem = np.append([freqMem],[readFreq[iLine/3]])
                    
                    tempMemIndex = np.append(tempMemIndex,len(tempMem))

                    if not np.count_nonzero(lineFreq == readFreq[iLine/3]): #If no choise has been done for this frequency
                        lineFreq = np.append(lineFreq,readFreq[iLine/3])
                        wave_length = ccgs / lineFreq[-1] * 1e7   #Wave length in nanometers
                        if wave_length < 10: bandNameInput = 'X-Rays'
                        elif (wave_length > 30000): 
                            if readFreqExp <= 11:  #Writing out numbers smaller than 1000 without power-of-ten
                                bandNameInputTemp = readoutFreq * 10 ** (readFreqExp-9)
                                if bandNameInputTemp == int(bandNameInputTemp):
                                    bandNameInputOut = int(readoutFreq * 10 ** (readFreqExp-9))
                                else:
                                    bandNameInputOut = readoutFreq * 10 ** (readFreqExp-9)
                                bandNameInput = '%s GHz'%bandNameInputOut
                            else:
                                bandNameInput = r'$%s \times 10^{%d}$ GHz'%(readoutFreq,readFreqExp-9)
                            
                        else:
                            bandNameInput = np.copy(bandNames[np.argmin(np.abs(bandWaveLengths - wave_length))])

                        if use_standard_colors:
                            linecolor = bandStandardColors[np.argmin(np.abs(standardFreqs-lineFreq[-1]))]
                        else:
                            linecolor = raw_input('Colour for %s (matplotlib standard, html hex or gray scale): '%(bandNameInput))
                        lineColorMem = np.append(lineColorMem,linecolor)


                        ### Scaling lightcurves

                        if (scaleLightcurves == 'y') or (scaleLightcurves == 'Y'):
                            if loadScaleData == '':
                                try: scaleLCinput = input("Scale factor for %s: "%bandNameInput)
                                except: scaleLCinput = input("Bad input! Scale factor for %s: "%bandNameInput)
                                saveScaleData = '%s%s %f\n'%(saveScaleData,lineFreq[-1],scaleLCinput)
                            else:
                                ### Finding the line in scaleData variable corresponging to the desired frequency
                                scaleDataIndex = np.argmin(np.abs(scaleData[:,0] - lineFreq[-1]))
                                if scaleData[scaleDataIndex,0] != lineFreq[-1]:
                                    new_scaleLCinput = raw_input('Data with scale factors has no input corresponding to the frequency %s. Type new scale factor  or leave blank if you want to use scale factor for the closest frequency (%s Hz): '%(lineFreq[-1],scaleData[scaleDataIndex,0]))
                                    if new_scaleLCinput != '':
                                        scaleLCinput = float(new_scaleLCinput)
                                    else:                                        
                                        scaleLCinput = scaleData[scaleDataIndex,1]
                                else:                                        
                                    scaleLCinput = scaleData[scaleDataIndex,1]

                            if scaleLCinput < 1: 
                                if (int(1./scaleLCinput) == (Decimal(str(1./scaleLCinput)).quantize(Decimal('0.01')))): #Removing numeric deviations
                                    bandName = np.append(bandName,'%s / %d'%(bandNameInput,int(1./scaleLCinput)))
                                elif ((int(1./scaleLCinput)+1) == (Decimal(str(1./scaleLCinput)).quantize(Decimal('0.01')))):
                                    bandName = np.append(bandName,'%s / %d'%(bandNameInput,int(1./scaleLCinput)+1))
                                else: bandName = np.append(bandName,'%s / %s'%(bandNameInput,1./scaleLCinput))
                            else: 
                                if float(scaleLCinput) == int(scaleLCinput):bandName = np.append(bandName,'%s x %d'%(bandNameInput,int(scaleLCinput)))
                                else: np.append(bandName,'%s x %d'%(bandNameInput,scaleLCinput))
                            scaleLC = np.append(scaleLC,np.copy(scaleLCinput))
                                
                                
                        else: 
                            scaleLC = np.append(scaleLC,1.)
                            bandName = np.append(bandName,bandNameInput)
                        legend[iLine/3] = bandName[-1]
                        
                    else: 
                        lineIndex = np.where(lineFreq==readFreq[iLine/3])
                        linecolor = lineColorMem[lineIndex][0]
                        lineColorMem = np.append(lineColorMem,linecolor)
                        scaleLC = np.append(scaleLC,scaleLC[lineIndex])


            ### Saving file with scale factors
            if (scaleLightcurves == 'y') or (scaleLightcurves == 'Y'):
                if (loadScaleData == ''):
                    saveScaleDataName = raw_input('Name of file to save scale factors in. Leave blank to not save: ')
                    if saveScaleDataName != '':
                        writeScaleData = open(saveScaleDataName,'w')
                        writeScaleData.write(saveScaleData)
                        writeScaleData.close()
                
            
            ### Printing legend with corresponding colors
            for iFreqLegend in range(len(lineFreq)):
                plt.plot(0,0,lineColorMem[iFreqLegend])
            plt.legend(legend,loc=input('Position of legend (1-4): '))


            for iPlot in range(len(freqMem)):
                plt.plot(tempMem[tempMemIndex[iPlot]:tempMemIndex[iPlot+1]]/scalePlotTime[daysOrSec],fluxMem[tempMemIndex[iPlot]:tempMemIndex[iPlot+1]] * scaleFluxAxis[fluxAxis] * scaleLC[iPlot],lineTypeMem[iPlot],color=lineColorMem[iPlot])

                ### Plotting uncertainty areas
                if plot_area:
                    try:
                        len(sigma_freq)
                        sigma_freq_is_scalar = False
                    except:
                        sigma_freq_is_scalar = True
                    if sigma_freq_is_scalar:
                        #which_sigma = np.where(sigma_freq == freqMem[iPlot])[0]
                        if (sigma_prob_input == 1) or (sigma_prob_input == 2):
                            plt.fill_between(sigmaTempGrid_68/scalePlotTime[daysOrSec] , highest_LC_68 * scaleFluxAxis[fluxAxis] * scaleLC[iPlot] , lowest_LC_68 * scaleFluxAxis[fluxAxis] * scaleLC[iPlot],color=color_68)
                        if (sigma_prob_input == 2):
                            plt.fill_between(sigmaTempGrid_95/scalePlotTime[daysOrSec] , highest_LC_95 * scaleFluxAxis[fluxAxis] * scaleLC[iPlot] , highest_LC_68 * scaleFluxAxis[fluxAxis] * scaleLC[iPlot],color=color_95)
                            plt.fill_between(sigmaTempGrid_95/scalePlotTime[daysOrSec] , lowest_LC_95 * scaleFluxAxis[fluxAxis] * scaleLC[iPlot] , lowest_LC_68 * scaleFluxAxis[fluxAxis] * scaleLC[iPlot],color=color_95)
                    else:
                        which_sigma = np.where(sigma_freq == freqMem[iPlot])[0]
                        if (sigma_prob_input == 1) or (sigma_prob_input == 2):
                            plt.fill_between(sigmaTempGrid_68[which_sigma][0]/scalePlotTime[daysOrSec] , highest_LC_68[which_sigma][0] * scaleFluxAxis[fluxAxis] * scaleLC[iPlot] , lowest_LC_68[which_sigma][0] * scaleFluxAxis[fluxAxis] * scaleLC[iPlot],color=color_68)
                        if (sigma_prob_input == 2):
                            plt.fill_between(sigmaTempGrid_95[which_sigma][0]/scalePlotTime[daysOrSec] , highest_LC_95[which_sigma][0] * scaleFluxAxis[fluxAxis] * scaleLC[iPlot] , highest_LC_68[which_sigma][0] * scaleFluxAxis[fluxAxis] * scaleLC[iPlot],color=color_95)
                            plt.fill_between(sigmaTempGrid_95[which_sigma][0]/scalePlotTime[daysOrSec] , lowest_LC_95[which_sigma][0] * scaleFluxAxis[fluxAxis] * scaleLC[iPlot] , lowest_LC_68[which_sigma][0] * scaleFluxAxis[fluxAxis] * scaleLC[iPlot],color=color_95)
            plot_area = False
                
    if printPlot:# ### Ordering input data. Useful when writing print-out and plotting data from print-out ###
        printFreqTemp = np.array(freq)
        printFreqOrd = np.zeros(iterationLength,dtype=int)
        
        for orderFreq in range(iterationLength):
            printFreqOrd[orderFreq] = np.argmin(printFreqTemp)
            printFreqTemp[printFreqOrd[orderFreq]] += 1e20

        printFreqOut = freq[printFreqOrd]
        tempGridOut = tempGrid[printFreqOrd]
        lightcurveOut = lightcurve[printFreqOrd]
        
        FdataPrint = np.copy(FdataInput[printFreqOrd])
        tdataPrint = np.copy(tdata[printFreqOrd])
        errorbarPrint = errorbarInput[printFreqOrd]

            
    for i in range(iterationLength):  #Looping over all frequencies


        ### Writing print-out

        if print_write: 
            printFreq, printFreqIte = printFreqOut[i], 0
            while printFreq >= 10:
                printFreq /= 10
                printFreqIte += 1

            
            readPrint += 'Band:%se%d\n%s\n'%(printFreq,printFreqIte,'\n'.join([' '.join(map(str,tempGridOut[i])),' '.join(map(str,lightcurveOut[i]))]))
            

        elif not print_write:   #If not read print-out
            if (plotDims == 'y') | (plotDims == 'Y'):
                plt.plot(freqGrid,FdataInput[:,0] * scaleFluxAxis[fluxAxis],'%so'%colourCycle[i%len(colourCycle)])
            
            if (not print_read) and  (not createMock): plt.plot(tempGrid[i]/scalePlotTime[daysOrSec],lightcurve[i] * scaleFluxAxis[fluxAxis],colourCycle[i%len(colourCycle)])



            ### Plotting 1-sigma area    ----

#            if plot_area:
#                for sigma_i in sigma_prob_range:
#                    if sigma_i == 0.68: plt.fill_between(sigmaTempGrid_68[i]/scalePlotTime[daysOrSec],highest_LC_68[i] * scaleFluxAxis[fluxAxis],lowest_LC_68[i] * scaleFluxAxis[fluxAxis],color=color_68)#,color=colourCycle[i%len(colourCycle)])
#                    elif sigma_i == 0.95: 
#                        plt.fill_between(sigmaTempGrid_95[i]/scalePlotTime[daysOrSec],highest_LC_95[i] * scaleFluxAxis[fluxAxis],highest_LC_68[i] * scaleFluxAxis[fluxAxis],color=color_95)
#                        plt.fill_between(sigmaTempGrid_95[i]/scalePlotTime[daysOrSec],lowest_LC_95[i] * scaleFluxAxis[fluxAxis],lowest_LC_68[i] * scaleFluxAxis[fluxAxis],color=color_95)#,color=colourCycle[i%len(colourCycle)])

                    


            if print_read: n_o_e_iterator = np.copy(numberOfEmpties[printFreqOrd[i]])
            else: n_o_e_iterator = np.copy(numberOfEmpties[i])
            for plotLC in range(n_o_e_iterator):
                    #Finding limits for plot
                    
                minFdataThis , maxFdataThis = FdataInput[i,np.argmin(FdataInput[i,:numberOfEmpties[i]])] , FdataInput[i,np.argmax(FdataInput[i,:numberOfEmpties[i]])]
                try: 
                    if minFdataThis < minFdata: minFdata = minFdataThis
                except: minFdata = minFdataThis
                try: 
                    if maxFdataThis > maxFdata: maxFdata = maxFdataThis
                except: maxFdata = maxFdataThis
                minTdataThis , maxTdataThis = tdata[i,np.argmin(tdata[i,:numberOfEmpties[i]])] , tdata[i,np.argmax(tdata[i,:numberOfEmpties[i]])]
                try: 
                    if minTdataThis < minTdata: minTdata = minTdataThis
                except: minTdata = minTdataThis
                try: 
                    if maxTdataThis > maxTdata: maxTdata = maxTdataThis
                except: maxTdata = maxTdataThis

                        
                if errorIsArr[i,plotLC]: #Plotting dot and errorbar. Else: plotting upper limit
                    if print_read:
                        plt.plot(tdata[printFreqOrd[i],plotLC]/scalePlotTime[daysOrSec],FdataInput[printFreqOrd[i],plotLC] * scaleFluxAxis[fluxAxis] * scaleLC[i],'%so'%lineColorMem[i])
                        plt.errorbar(tdata[printFreqOrd[i],plotLC]/scalePlotTime[daysOrSec],FdataInput[printFreqOrd[i],plotLC] * scaleFluxAxis[fluxAxis] * scaleLC[i],yerr=errorbarInput[i,plotLC],ecolor='%s'%lineColorMem[i])
#                        plt.plot(tdataPrint[i,plotLC]/scalePlotTime[daysOrSec],FdataPrint[i,plotLC] * scaleFluxAxis[fluxAxis] * scaleLC[i],'%so'%lineColorMem[i])
#                        plt.errorbar(tdataPrint[i,plotLC]/scalePlotTime[daysOrSec],FdataPrint[i,plotLC] * scaleFluxAxis[fluxAxis] * scaleLC[i],yerr=errorbarPrint[i,plotLC]*scaleLC[i],ecolor='%s'%lineColorMem[i])#,elinewidth=1)
                    else:
                        plt.plot(tdata[i,plotLC]/scalePlotTime[daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[fluxAxis],'%so'%colourCycle[i%len(colourCycle)])
                        plt.errorbar(tdata[i,plotLC]/scalePlotTime[daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[fluxAxis],yerr=errorbarInput[i,plotLC],ecolor='%s'%colourCycle[i%len(colourCycle)])#,elinewidth=1)
                else:
                    if print_read:
                        plt.errorbar(tdata[i,plotLC]/scalePlotTime[daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[fluxAxis]*0.75 * scaleLC[i],yerr=FdataPrint[i,plotLC] * scaleFluxAxis[fluxAxis]*0.25,lolims=True,ecolor='%s'%lineColorMem[i]) 
                    else:
                        plt.errorbar(tdata[i,plotLC]/scalePlotTime[daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[fluxAxis]*0.75,yerr=FdataInput[i,plotLC] * scaleFluxAxis[fluxAxis]*0.25,lolims=True,ecolor='%s'%colourCycle[i%len(colourCycle)]) 

                        
                #else:plt.errorbar(tdata[i,plotLC]/scalePlotTime[daysOrSec],FdataInput[i,plotLC] * scaleFluxAxis[fluxAxis]*0.75,yerr=FdataInput[i,plotLC] * scaleFluxAxis[fluxAxis]*0.25,lolims=True,ecolor='%s'%colourCycle[i%len(colourCycle)]) 

    if print_write:
        writePlotPrint = open(fileName,'w')
        writePlotPrint.write(readPrint)
        writePlotPrint.close()
    else:
        plt.xlim([0.5*minTdata/scalePlotTime[daysOrSec],2*maxTdata/scalePlotTime[daysOrSec]])
        plt.ylim([0.5*minFdata * scaleFluxAxis[fluxAxis],2*maxFdata * scaleFluxAxis[fluxAxis]])
                        
                
    #    elif mockDim == 'E': #option 'E' is obsolete
    #        plt.plot(freq[i,:numberOfEmpties[i]]/scalePlotTime,lightcurve[i,:numberOfEmpties[i]] * scaleFluxAxis,colourCycle[i%len(colourCycle)])
    #        if useData: plt.plot(freq[i,:numberOfEmpties[i]]/scalePlotTime,FdataInput[i,:numberOfEmpties[i]] * scaleFluxAxis,'%so'%colourCycle[i%len(colourCycle)])
    #    plt.loglog()
   # 

        if (plotDims == 'y') | (plotDims == 'Y'):  plt.xlabel('Frequency (Hz)') 
        else: plt.xlabel('Observing time (%s)'%('days'*(daysOrSec=='d') + 'hours'*(daysOrSec=='h') + 'minutes'*(daysOrSec=='m') + 'seconds'*(daysOrSec=='s')))
        plt.title(raw_input('Plot title: '))
        
#        if raw_input('Plot legend? (y/n)') == 'y':
#            legend = ['']*iterationLength
#            for iLegend in range(iterationLength):legend[iLegend] = 
#            plt.legend(legend)
        plt.ylabel('Flux (%s)'%fluxAxis)
        plt.loglog()
        plt.show()



    raise SystemExit(0)





#Produce surface plot of log-likelihood
if runOption == 'surf': 
    irange,jrange = np.arange(0,1,.05),np.arange(0,1,.05)
    xGrid,yGrid = np.zeros(len(irange)),np.zeros(len(jrange))
    loglike = np.zeros([len(irange),len(jrange)])

    #Loading mock observations
    try:loadMock = np.loadtxt('mockData.txt')
    except:
        if allowPrint: print "No previously saved mock observations (./mockData.txt) found!"
        raise SystemExit(0)
    tdata,Fdata = [],[]
    for imock in range(0,len(loadMock),2): 
        tdata.append(loadMock[imock])
        Fdata.append(loadMock[imock+1])
    errorbar = Fdata / 10. #Errorbar estimation for mock data


    for iSurfGrid in range(len(irange)):
        for jSurfGrid in range(len(jrange)):
            cube = [irange[iSurfGrid],jrange[jSurfGrid]]
            myPrior(cube,2,2)
            if iSurfGrid == 0: xGrid[jSurfGrid] = cube[1]
            if jSurfGrid == 0: yGrid[iSurfGrid] = cube[0]
            nl = len(freq)
            loglike[iSurfGrid][jSurfGrid] = logLikelihood(cube,n_params,n_params)
    xAxis,yAxis = np.meshgrid(xGrid,yGrid)
    tredfig = plt.figure()
    ax = Axes3D(tredfig)
    ax.plot_surface(xAxis,yAxis,np.log10(-loglike), rstride=1, cstride=1, cmap=plt.cm.hot)
    #ax.plot_surface(xAxis,yAxis,-loglike*2, rstride=1, cstride=1, cmap=plt.cm.hot)
    nameCount,j = 0,0
    titel = []
    for i in parametrar:
        if i: 
            titel.append(paramNames[j])
            nameCount += 1
        j += 1
    plt.xlabel(titel[1])
    plt.ylabel(titel[0])
    #plt.zlabel('Probability')
    plt.savefig('3Dplot.jpg')
    plt.show(ax)


    
    raise SystemExit(0)


elif runOption == 'marginal':
    plotMarginal(n_params)


############################################
#           Fitting routine                #
############################################

elif runOption == 'fit':
    createMock = False
    #Run MultiNest
    nl = len(freq)

    pymultinest.run(logLikelihood, myPrior, n_params, importance_nested_sampling = False, resume = True, verbose = True, sampling_efficiency = 'model', n_live_points = livePoints,evidence_tolerance=0.5)

    
    tidPost = np.array(time.gmtime())
    tidDiff = np.array(tidPost)-np.array(tidPre)
    tidMod = [24,60,60]
    if (tidDiff[1] != 0) & allowPrint: print "Crossing of month during simulation occured!"
    for tidCorr in range(5,2,-1): 
        if tidDiff[tidCorr] < 0: 
            tidDiff[tidCorr] += tidMod[tidCorr-3]
            tidDiff[tidCorr-1] -= 1
    if printCrucial: print "Fitting done! \n\n\nThat was exhausting, thank heaven it's done! Time elapsed: %dh %dm %ds. \n\nAwaiting orders!"%(tidDiff[3],tidDiff[4],tidDiff[5])


### Printing best fits ###

elif runOption=='print-stats':
     ### PRINTING STATS FROM chains/1-stats.dat ###
    path_to_chains = raw_input('Path to where to find chains directory. Use asterix analyse many paths: ')
    file_paths = glob.glob('%s/chains/1-stats.dat'%path_to_chains)
    column_title_input = np.zeros(len(file_paths),dtype='S256')
    outtext_first = '\\begin{tabular}{l'

    for i_column_titles in range(len(file_paths)):
        column_title_input[i_column_titles] = raw_input('Title for column of path %s: '%os.path.abspath(file_paths[i_column_titles]))
    chi2 = np.array([])
    chi2_sigma = np.array([])
    chi2_text = ''
    for i_file in range(len(file_paths)):
        readStats = open(file_paths[i_file],'r')
        statsIn = readStats.read().split('Dim No.')
        number_of_modes = int(statsIn[0].split('Total Modes Found:')[1].split()[0])
        print 'Posterior has %d mode%s'%(number_of_modes,'s'*(number_of_modes>1))
    #    gaussian_raw = statsIn[1].split('\n')
        number_of_parameters = np.sum(parametrar)
        
        gaussian = np.zeros([number_of_modes,number_of_parameters,2])
        best_fit = np.zeros([number_of_modes,number_of_parameters])
        stats_diff = np.zeros([number_of_modes,number_of_parameters,2])
        #if number_of_modes > 1: outtext = '\\begin{tabular}{l|l}\nMode 1 & \\\\\nParameter & Value\\\\\n'
        #else: outtext = '\\begin{tabular}{l|l}\nParameter & %s\\\\\n'%(' & '.join(column_title))
        round_factor = 2  #How many decimals
        
        for i_modes in range(number_of_modes):
            if i_file == 0 and i_modes == 0: column_title = '%s'%column_title_input[i_file]
            else: column_title += ' & %s'%column_title_input[i_file]

            if number_of_modes > 1:  column_title += ' (Mode %d)'%(i_modes+1)

            outtext_first += '|l'
#            if i_modes > 0: outtext += '%s Mode %d'%(column_titles[i_file],i_modes+1)
            if number_of_modes > 1:
                chi2 = np.append(chi2,-2*float(statsIn[0+i_modes*3].split('Strictly Local Log-Evidence')[1].split('+/-')[0].split()[-1]))
                chi2_sigma = np.append(chi2_sigma,2*float(statsIn[0+i_modes*3].split('Strictly Local Log-Evidence')[1].split('+/-')[1].split()[0]))
            else:
                chi2 = np.append(chi2,-2*float(statsIn[0].split('+/-')[0].split()[-1]))
                chi2_sigma = np.append(chi2_sigma,2*float(statsIn[0].split('+/-')[1].split()[0]))
            chi2_text += ' & $%s\\pm%s$'%(round_off(chi2[-1] / numberOfPoints,round_factor) , round_off(chi2_sigma[-1] / numberOfPoints,round_factor))
            for i_stats in range(number_of_parameters):
                gaussian[i_modes,i_stats] = np.array(map(float,statsIn[1+3*i_modes].split('\n')[i_stats+1].split()[1:]))
                best_fit[i_modes,i_stats] = float(statsIn[2+3*i_modes].split('\n')[i_stats+1].split()[1])
                if preferredScale[whereParam[i_stats]] == 'log':
                    stats_diff[i_modes,i_stats,0] = 10**(gaussian[i_modes,i_stats,0] + gaussian[i_modes,i_stats,1]) - 10**best_fit[i_modes,i_stats]
                    stats_diff[i_modes,i_stats,1] = 10**best_fit[i_modes,i_stats] - 10**(gaussian[i_modes,i_stats,0] - gaussian[i_modes,i_stats,1])
                    gaussian[i_modes,i_stats] = 10**gaussian[i_modes,i_stats]
                    best_fit[i_modes,i_stats] = 10**best_fit[i_modes,i_stats]
                else:
                    stats_diff[i_modes,i_stats,0] = gaussian[i_modes,i_stats,0] + gaussian[i_modes,i_stats,1] - best_fit[i_modes,i_stats]
                    stats_diff[i_modes,i_stats,1] = best_fit[i_modes,i_stats] - gaussian[i_modes,i_stats,0] + gaussian[i_modes,i_stats,1]
                if preferredScale[whereParam[i_stats]] == 'deg':
                    stats_diff[i_modes,i_stats,0] = 180 / np.pi * stats_diff[i_modes,i_stats,0]
                    stats_diff[i_modes,i_stats,1] = 180 / np.pi * stats_diff[i_modes,i_stats,1]
                    gaussian[i_modes,i_stats] = 180 / np.pi * gaussian[i_modes,i_stats]
                    best_fit[i_modes,i_stats] = 180 / np.pi * best_fit[i_modes,i_stats]

                if (stats_diff[i_modes,i_stats,0] < 0) or (stats_diff[i_modes,i_stats,1] < 0):
                    print 'best-fit value of parameter %s is outside of the standard deviation!'%paramNames[whereParam[i_stats]]
                    stats_diff[i_modes,i_stats] = np.abs(stats_diff[i_modes,i_stats])

                stats_fot = int(np.log10(best_fit[i_modes,i_stats])) #Factor of ten of best fit
                stats_fot -= (stats_fot < 0)  #Making sure we get the nearest factor of ten below the value
                if i_file == 0 and i_modes==0:
                    exec('outtext_line%d = \'%%s \'%%latexParamNamesLin[whereParam[i_stats]]'%(i_stats+i_modes*n_params))
#                if i_file == 0 and i_modes > 0:
#                    exec('outtext_line%d = \'%%s \'%%latexParamNamesLin[whereParam[i_stats]]'%(i_stats+i_modes*n_params))
                if (stats_fot >2) or (stats_fot < -2):
                    exec('outtext_line%d += \' & $%s^{+%s}_{-%s}\\\\times 10^{%d}$\''%(i_stats,round_off(best_fit[i_modes,i_stats]*10**(-stats_fot),round_factor) , round_off(stats_diff[i_modes,i_stats,0]*10**(-stats_fot),round_factor) , round_off(stats_diff[i_modes,i_stats,1]*10**(-stats_fot),round_factor) , stats_fot))
                else:
                    exec('outtext_line%d += \' & $%s^{+%s}_{-%s}$\''%(i_stats , round_off(best_fit[i_modes,i_stats],round_factor)  , round_off(stats_diff[i_modes,i_stats,0],round_factor) , round_off(stats_diff[i_modes,i_stats,1],round_factor)))
    outtext_lines = ''
    for join_lines in range(n_params):
        exec('outtext_lines += outtext_line%s'%join_lines)
        outtext_lines += '\\\\\n'
    outtext = '%s}\\\\\nParameters & %s\\\\\n%s$\\chi^2_{\\rm red}$ %s\\n\\end{tabular}'%(outtext_first,column_title,outtext_lines,chi2_text)
    output_file_name = raw_input('Output file name without suffix [output.tex]: ')
    if output_file_name == '': output_file_name = 'output'
    if os.path.isfile('%s.tex'%output_file_name):
        overwrite_input = raw_input('Overwrite file %s.tex? ([y]/n): '%output_file_name)
        if not ((overwrite_input == '') or (overwrite_input == 'y') or (overwrite_input == 'Y')):
            print 'Will not overwrite file %s.tex. Now exiting'%output_file_name
            raise SystemExit(0)
    write_outtext = open('%s.tex'%output_file_name,'w')
    write_outtext.write(outtext)
    write_outtext.close()
    

            
            


#############################################
#        Lightcurve production              #
#############################################


elif runOption == 'LC': #Producing light curves of choosen frequencies and constants from constants.txt
    if createMock and (runOption=='LC'): 
        FdataInput , tdata, errorbarInput, numberOfPoints  = [],[],[],0
    lightcurve_production(freq,numberOfMock,numberOfEmpties,FdataInput,tdata,errorbarInput)



 
#Prints passed time
