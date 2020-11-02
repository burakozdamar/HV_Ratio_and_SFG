#!/usr/bin/python
#coding=utf-8

import csv
import re
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
import linecache
from sys import argv
script, num, title = argv

# for f in os.listdir('.') :
	#if fnmatch.fnmatch(f, 'restricted_down_hbond_part*'):
	#	print(f)
	# if fnmatch.fnmatch(f, 'PEG2Water_part*'):
	# 	print(f)


a=1.3386
b=1.3286
#H-Horizontal, V-Vertical, N-Nagata, S-Simone
HS=6.75
VS=6.63
rS=1.02
def format(r):
    return "%.2f" % r
def ff(r):
    r = float(format(r))
    return r
def intra(f):
    intraL1=linecache.getline(f, 23)
    intraL1=intraL1[34:51]
    intraL1=float(intraL1)
    WL1=linecache.getline(f, 20)
    WL1=WL1[11:29]
    WL1=float(WL1)
    DintraL1= intraL1*WL1/(2*a*b)
    return DintraL1
def interL2(f):
    interL2=linecache.getline(f, 26)
    interL2=interL2[12:30]
    interL2=float(interL2)
    WL1=linecache.getline(f, 20)
    WL1=WL1[11:29]
    WL1=float(WL1)
    DinterL2= interL2*WL1/(a*b)
    return DinterL2
def interPEG(f):
    DinterPEG=linecache.getline(f, 3)
    DinterPEG=DinterPEG[39:58]
    DinterPEG=float(DinterPEG)
    return DinterPEG
def interref(f,f2):
    interR = interL2(f)+interPEG(f2)-VS
    interR = ff(interR)
    return interR
def intraref(f):
    intraR = intra(f) - HS
    intraR = ff(intraR)
    return intraR

def ratioR(f, f2):
    ratioR = (intra(f)/(interL2(f)+interPEG(f2))) - rS
    ratioR = ff(ratioR)
    return ratioR


#newFile = open("HVratioTOT.dat", "a+")


n=int(num)
Lratio=[]
for i in range(1,n+1):
	fHB = 'restricted_down_hbond_part_'+str(i)
	fPEG = 'PEG2Water_part_'+str(i)
	#newFile.write(str(i)+'  '+str(ratioR(fHB,fPEG))+'  '+str(ff(interPEG(fPEG)))+ '\n')
	Lratio.append(ratioR(fHB,fPEG))



# HN=7.01
# VN=6.71
# rN=1.05


def plotHBbars(H,V,name):
    label_list = ['part'+i for i in range(1,n+1)]
    x = range(1,n+1)

    intra = plt.bar([i+0.2 for i in x], H, width=0.4, alpha=0.8, color='red', label="Horizontal")
    inter = plt.bar([i+0.6 for i in x], V, width=0.4, alpha=0.8, color='blue', label="Vertical")
    plt.ylim(-4.0,4.0)
    plt.xlim(0.7,4)
    plt.xticks([index + 0.4 for index in x], label_list)
    plt.legend(prop={'size':7})
    plt.title(title+"-HV-ref-"+name, fontsize=10)
    plt.ylabel(r"HBs/nm2(SAM-refAir/Water)", fontsize=8)
    #Show values on each bar
    def add_values(data):
        for d in data:
            height = d.get_height()
            plt.text(d.get_x() + d.get_width() / 2,height+0.01, height, ha='center', va='bottom')
            d.set_edgecolor('white')

    #add_values(intra)
    #add_values(inter)
    

def plotRatio(r,name):
    label_list = ['part'+str(i) for i in range(1,n+1)]
    # label_list = ['part1', 'part2', 'part3']
    x = range(1,n+1)

    ratio = plt.bar([i+0.4 for i in x], r, width=0.4, alpha=0.8, color='black', label="ratio-ref-"+name)
    plt.ylim(-1.0,1.0)
    plt.xlim(0.7,n+1)
    plt.xticks([index + 0.4 for index in x], label_list)
    plt.title(title+"-ratio-ref-"+name)#, fontsize=10)
    plt.legend()#prop={'size':7})
    plt.ylabel(r"Ratio(SAM-refAir/Water)")#, fontsize=8)
    #plt.text(0.8, 3.8, "Ratio: {0}".format(r), size =8, color = "black", alpha = 1)
    #Show values on each bar
    def add_values(data):
        for d in data:
            height = d.get_height()
            plt.text(d.get_x() + d.get_width() / 2,height-0.1, height, ha='center', va='bottom')
            d.set_edgecolor('white')

    add_values(ratio)

#fig = plt.figure()
#ax1 = fig.add_subplot(2,2,1)
#plotHBbars(HrefN, VrefN, "Nagata")
#ax2 = fig.add_subplot(2,2,2)
#plotHBbars(HrefS, VrefS,  "Simone")
#ax3 = fig.add_subplot(2,2,3)
#plotRatio(ratioN ,"Nagata")
#ax4 = fig.add_subplot(2,2,4)


plotRatio(Lratio ,"Simone")
plt.tight_layout()
plt.savefig("RatioHVRefSimone.png")
plt.show()

