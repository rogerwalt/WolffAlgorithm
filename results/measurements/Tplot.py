import pylab as pl

import tdata as data

maxBinSingle = len(data.data[1]['results']['single']['usPerflip']['stderr']) - 1
maxBinWolff = len(data.data[1]['results']['wolff']['usPerflip']['stderr']) - 1

pltdataSingle = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['single']['usPerflip']['mean'], data.data[x]['results']['single']['usPerflip']['stderr'][maxBinSingle]) for x in range(len(data.data))]
pltdataSingle.sort(key=lambda x: (x[0], x[1]))

pltdataWolff = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['wolff']['usPerflip']['mean'], data.data[x]['results']['wolff']['usPerflip']['stderr'][maxBinWolff]) for x in range(len(data.data))]
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
    
pl.xlabel('T')
fig.suptitle('uS / spin flip', fontsize=14)
ax1.legend(bbox_to_anchor=(1.01, 0), loc='lower left', borderaxespad=0.)
pl.show()