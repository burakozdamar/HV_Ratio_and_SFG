import linecache
filename = 'nb_hbond'
L3=linecache.getline("./down_hbond", 4)
L3=L3[23:27]
print "Number of Hbonds in L3 : ", L3.rstrip()
#with open ('nb_hbond', 'w') as nb_hbond:
#    nb_hbond.write("L3hbonds = " L3.rstrip())

L2=linecache.getline("./down_hbond", 13)
L2=L2[23:27]
print "Number of Hbonds in L2 : ", L2.rstrip()
#with open ('nb_hbond', 'w') as nb_hbond:
#    nb_hbond.write("L2hbonds = " L2.rstrip())

L1=linecache.getline("./down_hbond", 22)
L1=L1[23:27]
#with open ('nb_hbond', 'w') as nb_hbond:
#     nb_hbond.write("L1hbonds = " L1.rstrip())
print "Number of Hbonds in L1 : ", L1.rstrip()

L0=linecache.getline("./down_hbond", 31)
L0=L0[23:27]
#with open ('nb_hbond', 'w') as nb_hbond:
#    nb_hbond.write("L0hbonds = " L0.rstrip())
print "Number of Hbonds in L0 : ", L0.rstrip()
