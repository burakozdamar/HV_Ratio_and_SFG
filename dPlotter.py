# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from decimal import Decimal

script, filename, line1, b1, b2= argv

line1=float(argv[2])
b1=float(argv[3])
b2=float(argv[4])

print(argv)

data = np.loadtxt(filename)

def local_min(ys):
        
    return [[i,np.arange(-4,25.2,0.2)[i],y] for i, y in enumerate(ys)
            if ((i == 0) or (ys[i - 1] >= y))
            and ((i == len(ys) - 1) or (y < ys[i+1]))]

print("Here are the minima!")
local_min_lst = local_min(data[:,1])
for i in range(len(local_min_lst)):
        print(local_min(data[:,1])[i]) 


first=line1
second=first+1.70
format(first,'1f')
format(second,'1f')
bound1, bound2 = -4, 12 #you can modify here to extend the boundaries of the plot

title2DN_lst=[first,second]
titleLayer_lst=[b1,b2]

distance_lst=list(np.arange(-400,2520,20)/100)

plt.plot(data[distance_lst.index(bound1):distance_lst.index(bound2),0],
         data[distance_lst.index(bound1):distance_lst.index(bound2),1],linewidth=2) 
plt.title("Restricted L1: {0}".format(title2DN_lst)+"  Layers Limits: {0}".format(titleLayer_lst))

plt.axvline(x=first, color='green', linewidth=1.2)
plt.axvline(x=second, color='green', linewidth=1.2)

plt.axvline(x=0.5, color='red', linewidth=1.2) 
plt.axvline(x=b1, color='red', linewidth=1.2) #uncomment this to draw a vertical line thru 1st minima
plt.axvline(x=b2, color='red', linewidth=1.2) #uncomment this to draw a vertical line thru 2nd minima

plt.ylim((0,2))
plt.xlabel(r'$r(\AA)$')
plt.ylabel(r'$\rho(r)/\rho(bulk)$')
plt.savefig(filename+".png",bbox_inches="tight",dpi=300)
plt.show()
