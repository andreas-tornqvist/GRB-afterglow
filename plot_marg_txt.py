def plot_marg_txt(param_names,whereParam):
    import numpy as np
    import os
    import glob
    from matplotlib import pyplot as plt

    print 'Available parameters: '
    for i in range(len(whereParam)):
        print '(%d): %s'%(i+1 , param_names[whereParam[i]])

    number_of_plots = input('Number of plots: ')
    try: number_of_columns = input('Number of subplot columns [%d]: '%int(np.sqrt(number_of_plots)))
    except: number_of_columns = int(np.sqrt(number_of_plots))

    pick_xaxis = np.zeros(number_of_plots,dtype=int)
    pick_yaxis = np.zeros(number_of_plots,dtype=int)

    while True:
        try: 
            exec('pick_xaxis[:] = [%s]'%raw_input('Enter the number to the x-axes, separated by commas: '))
            exec('pick_yaxis[:] = [%s]'%raw_input('Enter the number to the y-axes, separated by commas: '))
            for assign_pick in range(number_of_plots):
                pick_xaxis[assign_pick] = whereParam[pick_xaxis[assign_pick]]
                pick_yaxis[assign_pick] = whereParam[pick_yaxis[assign_pick]]
            break
        except:
            print 'Bad input! Try again'
        
    ### Finding how many sets of files there are for each parameter combination
    number_of_filesets = len(glob.glob('Figures/%s_%s_xaxis_*.txt'%(param_names[pick_xaxis[0]] , param_names[pick_yaxis[0]])))
    if number_of_filesets == 0: number_of_filesets = len(glob.glob('Figures/%s_%s_xaxis_*.txt'%(param_names[pick_yaxis[0]] , param_names[pick_xaxis[0]])))

    pick_fileset_input = raw_input('There are %d sets of files. Write which to load, separated by comma, or leave blank to load all: ')
    if pick_fileset_input == '': pick_fileset = range(number_of_filesets)
    else: exec('pick_fileset = [%s]'%pick_fileset_input)

    number_of_filesets = len(pick_fileset)

    this_figure , ax = plt.subplots( int(number_of_plots / number_of_columns) + ((int(number_of_plots / number_of_columns)*number_of_columns)<number_of_plots)    ,   number_of_columns)
    for i_plot in range(number_of_filesets): ### Looping over every set of files
        ### Reading in data from the file
        
        
    


        
