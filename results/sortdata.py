import data
import matplotlib.pyplot as plt

# sort data
data.data.sort(key=lambda tup: tup[0])


# plot cluster size
pltdata = [(data.data[x][0], data.data[x][1]['clusterSize']['mean'], data.data[x][1]['clusterSize']['stderr'][5]) for x in range(801)]
plt.errorbar([x for (x,y,z) in pltdata], [y for (x,y,z) in pltdata], [z for (x,y,z) in pltdata])

# plot energy
pltdata = [(data.data[x][0], data.data[x][1]['energy']['mean'], data.data[x][1]['energy']['stderr'][5]) for x in range(801)]

# plot magnetization
pltdata = [(data.data[x][0], data.data[x][1]['magnetization']['mean'], data.data[x][1]['magnetization']['stderr'][5]) for x in range(801)]
