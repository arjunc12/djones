import pylab
import os

# system calls for compiling and running code manually
os.system('gfortran funct.f -o f.x')
os.system('./f.x')

# load and plot data
data = pylab.loadtxt('plot.dat')

x = data[:,0]
yapprox = data[:,1]
yreal = pylab.e ** (-1 * x)

pylab.plot(x, yreal, c = 'r', label = 'e^-x')
pylab.plot(x, yapprox, c = 'b', linestyle = '--', label = 'forward euler')

pylab.show()