def bin_LC(read_name,print_name,n_points):
    import numpy as np
    indata = np.loadtxt(read_name)
    #n_points is number of points per bin
    a,b,c = binner(indata[:,0],n_points,'linear'),binner(indata[:,1],n_points,'linear'),binner(indata[:,2],n_points,'square')
    np.savetxt(print_name,np.transpose([a,b,c]))
    print 'File %s written'%print_name


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

read_name = raw_input('Read this file: ')
print_name = raw_input('Print this file: ')
n_points = input('Number of points per bin: ')
bin_LC(read_name,print_name,n_points)
