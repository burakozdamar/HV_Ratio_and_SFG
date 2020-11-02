"""This code checks if the water molecules are centered inside the sim. box"""
import numpy as np
import sys
with open(sys.argv[1]) as f:
    natom = int(f.readline())
    step = int(f.readline().split()[-1]) 
    lines = [next(f) for x in range(natom)]

z = [float(item.split()[-1]) for index, item in enumerate(lines) if index%3==0]

diff = abs(max(z)-abs(min(z)))

if diff < 1.:
    print("z_max: {:.2f}, z_min: {:.2f}, diff: {:.3f}".format(max(z), min(z), diff))
    print("The water molecules are reasonably centered in the simulation box. Kudos!")
else:
    print("Attention!")
    print("The water molecules do not seem to be centered in the simulation box. However, here is the z-coordinate distribution.")
    print("-"*16)
    print("z_max: {:.2f}, z_min: {:.2f}, diff: {:.3f}".format(max(z), min(z), diff))
    
    counts = np.histogram(z,bins=10)
   # print(counts)
    for i in range(len(counts[0])):
        #print("{:.2f},{:.2f} : {}".format (counts[1][i],counts[1][i+1],counts[0][i]))
        print("{:.2f},{:.2f} : {}".format (counts[1][-1-i],counts[1][-2-i],counts[0][-i-1]))
