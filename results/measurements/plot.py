import pylab as pl

import data_wolff_single_different_systemSizes_38e2c3474799a5e9cc375c760ab488327fdacf06 as data

xlbl='T'
maxBinSingle = len(data.data[1]['results']['single']['energy']['stderr']) - 1
maxBinWolff = len(data.data[1]['results']['wolff']['energy']['stderr']) - 1

# sort plotdata by systemsize and temperature
def plotComparison(measure):
    pltdataSingle = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['single'][measure]['mean'], data.data[x]['results']['single'][measure]['stderr'][maxBinSingle]) for x in range(len(data.data))]
    pltdataSingle.sort(key=lambda x: (x[0], x[1]))

    pltdataWolff = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['wolff'][measure]['mean'], data.data[x]['results']['wolff'][measure]['stderr'][maxBinWolff]) for x in range(len(data.data))]
    pltdataWolff.sort(key=lambda x: (x[0], x[1]))

    # initialize figure
    fig = pl.gcf()

    # initialize subplot
    ax1 = pl.subplot(211)

    for systemSize in sorted(list(set([a for (a,x,y,z) in pltdataSingle]))):
        pl.errorbar([x for (a,x,y,z) in pltdataSingle if a==systemSize],[y for (a,x,y,z) in pltdataSingle if a==systemSize],[z for (a,x,y,z) in pltdataSingle if a==systemSize], label=systemSize)

    pl.subplot(212, sharex = ax1)
    for systemSize in sorted(list(set([a for (a,x,y,z) in pltdataWolff]))):
        pl.errorbar([x for (a,x,y,z) in pltdataWolff if a==systemSize],[y for (a,x,y,z) in pltdataWolff if a==systemSize],[z for (a,x,y,z) in pltdataWolff if a==systemSize], label=systemSize)

    pl.xlabel(xlbl)
    fig.suptitle(title, fontsize=14)
    ax1.legend(bbox_to_anchor=(1.01, 0), loc='lower left', borderaxespad=0.)
    pl.show()

def plotAcceptanceRate():
    pltdataSingle = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['single']['acceptanceRate']['mean'], data.data[x]['results']['single']['acceptanceRate']['stderr'][maxBinSingle]) for x in range(len(data.data))]
    pltdataSingle.sort(key=lambda x: (x[0], x[1]))
    
    for systemSize in sorted(list(set([a for (a,x,y,z) in pltdataSingle]))):
                pl.errorbar([x for (a,x,y,z) in pltdataSingle if a==systemSize],[y for (a,x,y,z) in pltdataSingle if a==systemSize],[z for (a,x,y,z) in pltdataSingle if a==systemSize], label=systemSize)
    
    fig = pl.gcf()
    fig.suptitle(title,fontsize=14) 
    pl.xlabel(xlbl)
    pl.legend(loc='lower right')
    pl.show()


def plotClusterSize():
    pltdataSingle = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['wolff']['clusterSize']['mean'], data.data[x]['results']['wolff']['clusterSize']['stderr'][maxBinWolff]) for x in range(len(data.data))]
    pltdataSingle.sort(key=lambda x: (x[0], x[1]))

    for systemSize in sorted(list(set([a for (a,x,y,z) in pltdataSingle]))):
                pl.errorbar([x for (a,x,y,z) in pltdataSingle if a==systemSize],[y for (a,x,y,z) in pltdataSingle if a==systemSize],[z for (a,x,y,z) in pltdataSingle if a==systemSize], label=systemSize)

    fig = pl.gcf()
    fig.suptitle(title,fontsize=14)
    pl.xlabel(xlbl)
    pl.legend(loc='lower right')
    pl.show()
   

def plotClusterSizeRelative():
    pltdataSingle = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['wolff']['clusterSize']['mean'], data.data[x]['results']['wolff']['clusterSize']['stderr'][maxBinWolff]) for x in range(len(data.data))]
    pltdataSingle.sort(key=lambda x: (x[0], x[1]))

    for systemSize in sorted(list(set([a for (a,x,y,z) in pltdataSingle]))):
        pl.errorbar([x for (a,x,y,z) in pltdataSingle if a==systemSize],[y/(a**3) for (a,x,y,z) in pltdataSingle if a==systemSize],[z/(a**3) for (a,x,y,z) in pltdataSingle if a==systemSize], label=systemSize)

    fig = pl.gcf()
    fig.suptitle(title,fontsize=14)
    pl.xlabel(xlbl)
    pl.legend(loc='lower left')
    pl.show()

 
def plotComparisonAutoCorrelation(measure):
    pltdataSingle = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['single'][measure]['autocorr'][maxBinWolff]) for x in range(len(data.data))]
    pltdataSingle.sort(key=lambda x: (x[0], x[1]))

    pltdataWolff = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['wolff'][measure]['autocorr'][maxBinWolff]) for x in range(len(data.data))]
    pltdataWolff.sort(key=lambda x: (x[0], x[1]))

    # initialize figure
    fig = pl.gcf()

    # initialize subplout
    ax1 = pl.subplot(211)

    for systemSize in sorted(list(set([a for (a,x,y) in pltdataSingle]))):
        pl.plot([x for (a,x,y) in pltdataSingle if a==systemSize],[y for (a,x,y) in pltdataSingle if a==systemSize], label=systemSize)

    pl.subplot(212, sharex = ax1)
    for systemSize in sorted(list(set([a for (a,x,y) in pltdataWolff]))):
        pl.plot([x for (a,x,y) in pltdataWolff if a==systemSize],[y for (a,x,y) in pltdataWolff if a==systemSize], label=systemSize)

    pl.xlabel(xlbl)
    fig.suptitle('Autocorrelation: '+title, fontsize=14)
    ax1.legend(bbox_to_anchor=(1.01, 0), loc='lower left', borderaxespad=0.)
    pl.show()

# plot cluster size
title='Relative Cluster size'
plotClusterSizeRelative()

# plot autocorrelations
title='Energy'
plotComparisonAutoCorrelation('energy')

title='Magnetization'
plotComparisonAutoCorrelation('magnetization')

# plot energy
title='Energy'
plotComparison('energy')

# plot magnetization
title='Magnetization'
plotComparison('magnetization')

# plot acceptance rate
title='Acceptance rate'
plotAcceptanceRate()

# plot cluster size
title='Cluster size'
plotClusterSize()
