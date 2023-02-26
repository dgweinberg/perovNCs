
import numpy as np
import os
import sys

filename = sys.argv[1]
if not os.path.isfile(filename):
   print ('File name {} does not exist.'.format(filename))
   sys.exit()

output=open('btraj.xyz', 'w')

print ('')

dt = 1.0         # timestep
nframe = 0     # of frames (snapshots)

with open(filename) as inp:
   line = inp.readline()
   while len(line) != 0:
       # converting lammps dump file to VMD readable xyz format
       nframe += 1
       nch = 0
       tstep = inp.readline()
       line = inp.readline()
       if nframe % 100 == 1: print ('%i th frame' %nframe)
       if nframe == 1:
          tstart = int(tstep)*dt
       tend = int(tstep)*dt
       natom = inp.readline()
       line = inp.readline()
       xlo, xhi = map(float, inp.readline().split())
       ylo, yhi = map(float, inp.readline().split())
       zlo, zhi = map(float, inp.readline().split())
       # ----------------------------
       ntot = int(natom)
       line = inp.readline()
       atype = np.zeros(ntot)
       xu = np.zeros(ntot)
       yu = np.zeros(ntot)
       zu = np.zeros(ntot)
       q = np.zeros(ntot)
       nout = 0	# of outside atoms
       npb = 0		# of pb atoms
       ncs = 0
       ni = 0
       for i in range(ntot):
           atype[i], xu[i], yu[i], zu[i], q[i] = map(float, inp.readline().split())
           if atype[i] == 3: ncs += 1
           if atype[i] == 1: npb += 1
           if atype[i] == 2: ni += 1
           #if nframe == 1 and i < 10: print (i+1, atype[i], xu[i], yu[i], zu[i], q[i])

       # writing vmd file
       nstep = 1
       if nframe % nstep == 0: output.write('%i \n' %ntot)
       if nframe % nstep == 0: output.write('Atoms \n')
       for i in range(ntot):
           nch += 1
           if int(atype[i]) > 5: print ('check atom type')
           if int(atype[i]) == 1:
              if nframe % nstep == 0: output.write('Pb ')
           elif int(atype[i]) == 2:
              if nframe % nstep == 0: output.write('I ')
           elif int(atype[i]) in [3]:
              if nframe % nstep == 0: output.write('Cs ')
           if nframe % nstep == 0: output.write('  %lf' %xu[i])
           if nframe % nstep == 0: output.write('  %lf' %yu[i])
           if nframe % nstep == 0: output.write('  %lf \n' %zu[i])
           #if nframe == 1 and i < 10: print (i+1, atype[i], xu[i], yu[i], zu[i])
       line = inp.readline()
       if nch != len(atype)+0:
          print ('check array length again')
          print (nch, len(atype)+0)

print ('')
print ('nstep = %i ' %nstep)
print ('number of frames = %i ' %nframe)
print ('')
#
conv = np.power(1.0, 0.0)
print ('')
print ('************************ total simulation time *************************')
print ('tstart = %lf' %(tstart/1000.), 'ps', ' & tend = %lf' %(tend/1000.), 'ps')
print ('')



