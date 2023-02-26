
import numpy as np
import os
import sys

print ('')

dt = 1.0         # timestep
nframe = 0     # of frames (snapshots)
latc = []
latcx = []
latcy = []
latcz = []
latcd = []
angle = []
bond = []

# average LC (all three directions)
latc.append(6.442624)
latc.append(6.441743)
latc.append(6.441266)
latc.append(6.441201)
latc.append(6.440589)


# LC diff. from all coordinates
latcd.append(0.039902)
latcd.append(0.024437)
latcd.append(0.030743)
latcd.append(0.026743)
latcd.append(0.021209)

# LC diff. from pb distance
#latcd.append(0.100658)
#latcd.append(0.076745)
#latcd.append(0.077656)
#latcd.append(0.098682)
#latcd.append(0.100600)

angle.append(156.301204)
angle.append(156.051967)
angle.append(156.034984)
angle.append(155.825617)
angle.append(156.290328)

bond.append(3.204322)
bond.append(3.204824)
bond.append(3.205208)
bond.append(3.205298)
bond.append(3.205177)


print ('')
avglc = np.mean(latc)
stdlc = np.std(latc)
print ('avg. 1/2 LC = %lf' %(avglc*0.5))
print ('avg, LC = %lf' %(avglc))
print ('std. LC = %lf' %(stdlc))
print ('')
avglcd = np.mean(latcd)
stdlcd = np.std(latcd)
print ('average LC diff. = %lf' %avglcd)
print ('std of LC diff. = %lf' %stdlcd)
print ('')
avgag = np.mean(angle)
stdag = np.std(angle)
print ('avg angle = %lf' %avgag)
print ('std angle = %lf' %stdag)
print ('')
avgb = np.mean(bond)
stdb = np.std(bond)
print ('avg bond = %lf' %avgb)
print ('std bond = %lf' %stdb)
print ('')
#print ('************************ total simulation time *************************')
#print ('tstart = %lf' %(tstart/1000.), 'ps', ' & tend = %lf' %(tend/1000.), 'ps')
#print ('')



