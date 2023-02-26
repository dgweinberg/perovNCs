
import numpy as np
import os
import sys

filename = sys.argv[1]
if not os.path.isfile(filename):
   print ('File name {} does not exist.'.format(filename))
   sys.exit()

print ('')

dt = 1.0         # timestep
nframe = 0     # of frames (snapshots)
latc = []
latcx = []
latcy = []
latcz = []
latdiff = []

with open(filename) as inp:
   line = inp.readline()
   while len(line) != 0:
       # converting lammps dump file to VMD readable xyz format
       nframe += 1
       nch = 0
       tstep = inp.readline()
       line = inp.readline()
       #if nframe % 100 == 1: print ('%i th frame' %nframe)
       if nframe == 1:
          tstart = int(tstep)*dt
       tend = int(tstep)*dt
       natom = inp.readline()
       line = inp.readline()
       xlo, xhi = map(float, inp.readline().split())
       ylo, yhi = map(float, inp.readline().split())
       zlo, zhi = map(float, inp.readline().split())
       dlx = xhi - xlo
       dly = yhi - ylo
       dlz = zhi - zlo
       mindl = min(dlx, dly, dlz)
       maxdl = max(dlx, dly, dlz)
       if nframe == 1: print ('box length diff. = %lf' %(maxdl - mindl))
       # ----------------------------
       ntot = int(natom)
       line = inp.readline()
       atype = np.zeros(ntot)
       xu = np.zeros(ntot)
       yu = np.zeros(ntot)
       zu = np.zeros(ntot)
       q = np.zeros(ntot)
       xpb = []
       ypb = []
       zpb = []
       nout = 0	# of outside atoms
       npb = 0		# of pb atoms
       ncs = 0
       ni = 0
       for i in range(ntot):
           atype[i], xu[i], yu[i], zu[i], q[i] = map(float, inp.readline().split())
           if atype[i] == 3: ncs += 1
           if atype[i] == 1: npb += 1
           if atype[i] == 2: ni += 1
           if atype[i] == 1: xpb.append(xu[i])
           if atype[i] == 1: ypb.append(yu[i])
           if atype[i] == 1: zpb.append(zu[i])

       if nframe == 1: print ('')
       if nframe == 1: print ('frame # = %i' %nframe)
       ncube = int(np.power(npb+0.1, 1./3.))
       blx = (max(xu) - min(xu))/float(ncube-0.5)
       bly = (max(yu) - min(yu))/float(ncube-0.5)
       blz = (max(zu) - min(zu))/float(ncube-0.5)
       minbl = min(blx, bly, blz)
       maxbl = max(blx, bly, blz)
       if nframe == 1: print ('LC diff from all coord = %lf' %(maxbl - minbl))
       if nframe == 1: print ('cube # = %i' %ncube)
       latpbx = (max(xpb)-min(xpb))/float(ncube-1.)
       latpby = (max(ypb)-min(ypb))/float(ncube-1.)
       latpbz = (max(zpb)-min(zpb))/float(ncube-1.)
       minlc = min(latpbx, latpby, latpbz)
       maxlc = max(latpbx, latpby, latpbz)
       #latdiff.append(maxlc-minlc)
       if nframe == 1: print ('lattice constant (x) from Pb = %lf ' %latpbx)
       if nframe == 1: print ('lattice constant (y) from Pb = %lf ' %latpby)
       if nframe == 1: print ('lattice constant (z) from Pb = %lf ' %latpbz)
       if nframe == 1: print ('lattice constant (diff) from Pb = %lf ' %(maxlc-minlc))
       if nframe == 1: print ('')
       latc.append(blx)
       latc.append(bly)
       latc.append(blz)
       latcx.append(blx)
       latcy.append(bly)
       latcz.append(blz)
       latdiff.append(maxbl-minbl)

       #latc.append(latpbx)
       #latc.append(latpby)
       #latc.append(latpbz)      
       #latcx.append(latpbx)
       #latcy.append(latpby)
       #latcz.append(latpbz)
       line = inp.readline()

print ('')
avglc = np.mean(latc)
stdlc = np.std(latc)
avglcx = np.mean(latcx)
avglcy = np.mean(latcy)
avglcz = np.mean(latcz)
avgdiff = np.mean(latdiff)
print ('1/2 average LC = %lf' %(avglc*0.5))
print ('average LC = %lf' %avglc)
print ('std of LC = %lf' %stdlc)
print ('')
print ('avg. (x) LC = %lf' %avglcx)
print ('avg. (y) LC = %lf' %avglcy)
print ('avg. (z) LC = %lf' %avglcz)
print ('')
print ('avg. LC diff = %lf' %avgdiff)
#print ('number of frames = %i ' %nframe)
#print ('')
#
conv = np.power(1.0, 0.0)
print ('')
#print ('************************ total simulation time *************************')
#print ('tstart = %lf' %(tstart/1000.), 'ps', ' & tend = %lf' %(tend/1000.), 'ps')
#print ('')



