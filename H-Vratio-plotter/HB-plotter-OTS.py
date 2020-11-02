#!/usr/bin/python
#coding=utf-8

import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import linecache
from sys import argv
script, title = argv 

fana1 = "restricted_down_hbond_part_1" 
fana2 = "restricted_down_hbond_part_2" 
fana3 = "restricted_down_hbond_part_3" 


a=1.3386
b=1.3286


def intra(f):
    intraL1=linecache.getline(f, 23)
    intraL1=intraL1[34:51]
    intraL1=float(intraL1)
    WL1=linecache.getline(f, 20)
    WL1=WL1[11:29]
    WL1=float(WL1)
    DintraL1= intraL1*WL1/(2*a*b)
    return DintraL1
def inter(f):
    interL2=linecache.getline(f, 26)
    interL2=interL2[12:30]
    interL2=float(interL2)
    WL1=linecache.getline(f, 20)
    WL1=WL1[11:29]
    WL1=float(WL1)
    DinterL2= interL2*WL1/(a*b)
    return DinterL2
def format(r):
    return "%.2f" % r
def ff(r):
    r = float(format(r)) 
    return r
def intraref(f, ref):
    intraR = intra(f) - ref
    return intraR
def interref(f, ref):
    interR = inter(f)- ref
    interR = float(interR)
    return interR
def ratioR(f, ref):
    ratioR = (intra(f)/inter(f)) - ref     
    return ratioR
#H-Horizontal, V-Vertical, N-Nagata, S-Simone
HN=7.01
HS=6.75
VN=6.71
VS=6.63
rN=1.05
rS=1.02

HrefN = [intraref(fana1, HN), intraref(fana2, HN), intraref(fana3, HN)]
VrefN = [interref(fana1, VN), interref(fana2, VN), interref(fana3, VN)]
ratioN = [ratioR(fana1, rN), ratioR(fana2, rN), ratioR(fana3, rN)]

HrefS = [intraref(fana1, HS), intraref(fana2, HS), intraref(fana3, HS)]
VrefS = [interref(fana1, VS), interref(fana2, VS), interref(fana3, VS)]
ratioS = [ratioR(fana1, rS), ratioR(fana2, rS), ratioR(fana3, rS)]

ratioN = [ff(i) for i in ratioN]
HrefN = [ff(i) for i in HrefN]
VrefN = [ff(i) for i in VrefN]

ratioS = [ff(i) for i in ratioS]
HrefS = [ff(i) for i in HrefS]
VrefS = [ff(i) for i in VrefS]


def plotHBbars(H,V,name):
    label_list = ['part1', 'part2', 'part3']
    x = range(1,4)

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
    label_list = ['part1', 'part2', 'part3']
    x = range(1,4)

    ratio = plt.bar([i+0.4 for i in x], r, width=0.4, alpha=0.8, color='black', label="ratio-ref-"+name)
    plt.ylim(-0.8,0.5)
    plt.xlim(0.7,4)
    plt.xticks([index + 0.4 for index in x], label_list)
    plt.title(title+"-ratio-ref-"+name)
    #plt.title(title+"-ratio-ref-"+name, fontsize=10)
    plt.legend()
    #plt.legend(prop={'size':7})
    plt.ylabel(r"Ratio(SAM-refAir/Water)")
    #plt.ylabel(r"Ratio(SAM-refAir/Water)", fontsize=8)
    #plt.text(0.8, 3.8, "Ratio: {0}".format(r), size =8, color = "black", alpha = 1)
    #Show values on each bar
    def add_values(data):
        for d in data:
            height = d.get_height()
            plt.text(d.get_x() + d.get_width() / 2, height-0.1, height, ha='center', va='bottom')
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
plotRatio(ratioS ,"Simone")
plt.tight_layout()
plt.savefig('RatioHVRefSimone.png')
plt.show()


