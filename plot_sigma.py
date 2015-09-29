#This function reads in the posterior points, finds the 68% range and feeds it to the lightcurve generator, in order to plot the 1-sigma field on the lightcurve plot
def binner(x_array,y_array,bins):
    import numpy as np
    x_output , y_output = np.zeros(bins) , np.zeros(bins)
    bin_space = np.linspace(x_array[0],x_array[-1],bins+1)

    for i in range(bins):
        lower_lim = np.argmin(np.abs(x_array - bin_space[i]))
        upper_lim = np.argmin(np.abs(x_array - bin_space[i+1]))
        if lower_lim == upper_lim:
            x_output[i] = x_array[lower_lim]
            y_output[i] = 0
            continue
        else:
            x_output[i] = np.sum(x_array[lower_lim:upper_lim]) / (upper_lim - lower_lim)   #summing all elements  lower_lim <= i < upper_lim and normalizing
            y_output[i] = np.sum(y_array[lower_lim:upper_lim])    #Summing all elements, without normalizing
    
    return x_output , y_output , bin_space


def get_one_sigma(n_params,sigma_range):

    from fitterFunction import modelFunc
    from options import userOptions
    import numpy as np
    from matplotlib import pyplot as plt
    from numpy.core.defchararray import split as splt

#Read in chains/1-post_separate.dat

    file_name = 'chains/1-post_separate.dat'

#Finding number of modes

    plot_area = True

    try:
        read_text = open(file_name,'r')
        dat_text = read_text.read().split('\n')
        read_text.close()
        raw_text = np.array(dat_text)
    except:
        print 'No such file %s. Now exiting'%file_name
        raise SystemExit(0)

    empties = np.where(raw_text=='')[0]
    modes = 0

    for i in empties:
        try:
            if raw_text[i+1] == '': modes += 1
        except:continue
    print '%s has %d nodes'%(file_name,modes)

    in_data = np.loadtxt(file_name)
    length = len(in_data)

    tot_prob = 0.
    
    which_points = np.array([],dtype=int)
    order = np.argsort(-in_data[:,0]) #Sorted on probability, highest first
    sorted_prob = in_data[order,0]

    for i in range(length):
        tot_prob += sorted_prob[i]
        which_points = np.append(which_points,order[i])
        if tot_prob >= sigma_range: break
    
    ### Finding difference in total chi2
    
        
    highest_chi2 , lowest_chi2 = np.max(in_data[which_points,1]) , np.min(in_data[which_points,1])
    print 'Highest chi2 in distribution: %s\nLowest chi2 in distribution: %s\nchi2 difference: %s'%(highest_chi2,lowest_chi2,highest_chi2-lowest_chi2)
    
    return in_data[which_points]




