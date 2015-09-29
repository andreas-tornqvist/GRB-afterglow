import numpy as np

readStats = open('chains/1-stats.dat','r')
statsIn = readStats.read().split('Sigma')[1]
stats = statsIn.split('\n')[1:]
#print stats
#print stats[1].split()

while stats[0].split()[0] != '1':
#    print float(stats[0].split()[0])
    stats = stats[1:]

n=0
string = map(float,stats[n].split())
nu_mean = np.array([])
sigma = np.array([])
oneMore = True
while oneMore:
    nu_mean = np.append(nu_mean,string[1])
    sigma = np.append(sigma,string[2])
    n+=1
    string = map(float,stats[n].split())
    #print string
    #print n
    try:
        if string[0] == (n+1): oneMore = True
        else: oneMore = False
     #   print string[0]
    except: oneMore = False
print ' '.join(map(str,nu_mean))
print ' '.join(map(str,sigma))
