import pylab as pl

import tdata as data

pltdataSingle = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['single']['nsPerflip']) for x in range(len(data.data))]
pltdataSingle.sort(key=lambda x: (x[0], x[1]))

pltdataWolff = [(data.data[x]['systemSize'], data.data[x]['temperature'], data.data[x]['results']['wolff']['nsPerflip']) for x in range(len(data.data))]
pltdataWolff.sort(key=lambda x: (x[0], x[1]))

# initialize figure
fig = pl.gcf()

# initialize subplot
ax1 = pl.subplot(211)

for systemSize in sorted(list(set([a for (a,x,y) in pltdataSingle]))):
    pl.plot([x for (a,x,y) in pltdataSingle if a==systemSize],[y for (a,x,y) in pltdataSingle if a==systemSize], label=systemSize)

pl.subplot(212, sharex = ax1)

for systemSize in sorted(list(set([a for (a,x,y) in pltdataWolff]))):
    pl.errorbar([x for (a,x,y) in pltdataWolff if a==systemSize],[y for (a,x,y) in pltdataWolff if a==systemSize], label=systemSize)
    
pl.xlabel('T')
fig.suptitle('ns / spin flip', fontsize=14)
ax1.legend(bbox_to_anchor=(1.01, 0), loc='lower left', borderaxespad=0.)
pl.show()