#Palaiseau ECLA, 2020.07.14
#Wanlin
#SFG_plotter in python to replace xmGrace XD

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re
from sys import argv

script, filename = argv
disData = np.loadtxt(filename)
disK = disData[:, 1]
disCl = disData[:, 2]
disKCl = disData[:, 3]

def plotDis(disWhich, ionType, colorWhich):
	x=disData[:, 0]
	plt.plot(x,disWhich,label=ionType,color = colorWhich)
	plt.legend()
	plt.title(filename[0:len(filename)-4])
	#plt.xlim(0,65)
	plt.ylim(0, 15)
	plt.xlabel(r"t(ps)")
	plt.ylabel(r"$Distance(\AA)$")


plotDis(disK, "K+", "red")
plotDis(disCl, "Cl-", "blue")
plotDis(disKCl, "K-Cl", "green")
plt.savefig(filename[0:len(filename)-4]+".png")
plt.show()



