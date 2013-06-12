import data
import matplotlib.pyplot as plt

# sort data
d.data.sort(key=lambda tup: (tup['temperature'],tup['systemSize']))

# plot cluster size
pltdata = [(d.data[x]['temperature'], d.data[x]['results']['wolff']['energy']['mean']) for x in range(len(d.data))]
plt.errorbar([x for (x,y,z) in pltdata], [y for (x,y,z) in pltdata], [z for (x,y,z) in pltdata])

# plot energy
pltdata = [(data.data[x][0], data.data[x][1]['energy']['mean'], data.data[x][1]['energy']['stderr'][5]) for x in range(801)]

# plot magnetization
pltdata = [(data.data[x][0], data.data[x][1]['magnetization']['mean'], data.data[x][1]['magnetization']['stderr'][5]) for x in range(801)]
