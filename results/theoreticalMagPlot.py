import numpy
import matplotlib.pyplot as plt

x=numpy.arange(0,1,0.001)
y=numpy.power(1-x,0.326)

x=numpy.append(x,1)
y=numpy.append(y,0)

x=numpy.append(x,2)
y=numpy.append(y,0)

plt.ylim((-0.5, 1.5))

fig=plt.gcf()
fig.suptitle('Magnetization',fontsize=14)
plt.xlabel('T / Tc')
plt.plot(x,y);plt.show()
