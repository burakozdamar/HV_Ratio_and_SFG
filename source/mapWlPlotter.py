#2020.06.17
#Palaiseau ECLA
#Wanlin Chen 

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from sys import argv
from matplotlib import axes

script, filename= argv

data = np.loadtxt(filename)
a, b = 13.386, 13.286
nlx, nly = 13, 13
dlx = a/nlx
dly = b/nly
ai = -a/2
bi = b/2

# row 1 column 1 integration of 4 parts
R1C1 = np.sum(data[:6, :6])
R1C2 = np.sum(data[:6, 6:])
R2C1 = np.sum(data[6:, :6])
R2C2 = np.sum(data[6:, 6:]) 

#print(R1C1,R1C2, R2C1, R2C2)

def format(r):
    return "%.1f" % r 
def ff(r):
    float(format(r))

def draw():
    xLabel = [format(ai+i*dlx) for i in range(13)]
    yLabel = [format(bi-i*dly) for i in range(13)]
    #plt.xlim(ai, ai+a)
    #plt.ylim(bi, bi+b)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_yticks(range(13))
    ax.set_yticklabels(yLabel)
    ax.set_xticks(range(13))
    ax.set_xticklabels(xLabel)
    #plot heat plot
    #map = plt.imshow(data, cmap='jet')#GnBu')
    map = plt.imshow(data, interpolation='spline16', cmap='jet')
    #map = plt.imshow(data, vmin=0.0000, vmax = 0.4000, interpolation='spline16', cmap='jet')
    #colorbar on the right
    plt.colorbar(map)
    #plt.text(3,3,format(R1C1))
    #plt.title(filename[4:len(filename)-4])
    plt.savefig(filename[0:(len(filename)-3)]+"png")
#    plt.show()

draw()

