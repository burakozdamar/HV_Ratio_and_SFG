#Palaiseau ECLA, 2020.07.14
#Wanlin
#SFG_plotter in python to replace xmGrace XD

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re
from sys import argv

script, filename1,outname = argv

#directory = "."
N=50
n=np.ones(N)
weights=n/N

	#filename_re = Lx+'_part_\d'
	#for filename in os.listdir(directory):
		#if re.search(filename_re, filename):
Data_part1 = np.loadtxt(filename1)


def plotSFG(Data_Spectrum, colors):
	#Data_Spectrum=np.loadtxt('L1_part7to12.dat')
	N=50
	n=np.ones(N)
	weights=n/N
	#x = np.convolve(weights,Data_Spectrum[1:4000,0])[N-1:-N+1]
	y = np.convolve(weights,-Data_Spectrum[1:4000,2])[N-1:-N+1]
	#y = -Data_Spectrum[:,2]
	x=np.arange(len(y))
	plt.plot(x,y,color = colors)
	plt.legend()
	plt.xlim(2800,4000)
	plt.ylim(-10, 10)
	plt.xlabel(r"$Frequency(cm^-1)$")
	plt.ylabel(r"$\chi_{ssp}^2$"+r"$(10^{-22}m^2V^{-1})$")
	#plt.show()

plotSFG(Data_part1,'black')
#plotSFG(Data_part2,'blue',name2)
#plotSFG(Data_part3,'red',name3)
yy = [0,0,0]
xx =[2800,3500,4000]
plt.plot(xx,yy,color='grey')
plt.savefig(outname)
plt.show()


