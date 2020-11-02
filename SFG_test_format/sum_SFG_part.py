#Evry,2020.07.10,Wanlin Chen

import os
import re
import numpy as np
from sys import argv

script, Ni, Nf = argv 

#directory = "/Volumes/data-wanlin/ANA_SAM/OTS8_traj1/ana_tot/SFG_test/total"

#filename_re = 'ssp_dL1_\d'
#tot_array = np.zeros((4000,3))
#nbfile = 0
#for filename in os.listdir(directory):
#    if re.search(filename_re, filename):
#	    print(filename)
#	    nbfile += 1  
#	    temp = np.loadtxt(filename)
#	    tot_array += temp
#
#tot_array /= nbfile
#    np.savetxt('L1_tot_test.dat', tot_array)
#
def partAverage(a, b, layer):
    nb = 0 
    part_array = np.zeros((4000,3))
    for i in range(a,b+1):
        filepart = 'ssp_'+layer+'_'+str(i)
        data= np.loadtxt(filepart)
        nb += 1
        part_array += data 
        print(filepart)
    part_array /= nb 
    np.savetxt(layer+'_part'+str(a)+'to'+str(b)+'.dat',part_array)
Ni = int(Ni)
Nf = int(Nf)
#partAverage(Ni, Nf, 'L0')
partAverage(Ni, Nf, 'dL1')
partAverage(Ni, Nf, 'dL2')
partAverage(Ni, Nf, 'dL3')
partAverage(Ni, Nf, 'down')
partAverage(Ni, Nf, 'up')
