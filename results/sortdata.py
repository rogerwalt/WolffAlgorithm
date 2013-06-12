import data
import matplotlib.pyplot as plt

# sort data
data.data.sort(key=lambda tup: tup[0])


# plot cluster size
pltdata = [(data.data[x][0], data.data[x][1]['clusterSize']['mean'], data.data[x][1]['clusterSize']['stderr'][5]) for x in range(801)]
plt.plot([x for (x,y,z) in pltdata], [y for (x,y,z) in pltdata])
