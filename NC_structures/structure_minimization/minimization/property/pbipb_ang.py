
import numpy as np
import os
import sys
import math

filename = sys.argv[1]
if not os.path.isfile(filename):
   print ('File name {} does not exist.'.format(filename))
   sys.exit()

output1=open('avg_pbipb.dat', 'w')
output2=open('sur_pbipb.dat', 'w')

dt = 1.0
nframe = 0     # of frames (snapshots)
print ('')
print ('# compute distribution for PbI distance #')
print ('')

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
       # ----------------------------
       ntot = int(natom)
       line = inp.readline()
       atype = np.zeros(ntot)
       xu = np.zeros(ntot)
       yu = np.zeros(ntot)
       zu = np.zeros(ntot)
       q = np.zeros(ntot)
       nout = 0   # of outside atoms
       npb = 0    # of pb atoms
       ncs = 0
       ni = 0
       qtot = 0.0
       xpb = []
       indpb = []
       for i in range(ntot):
           atype[i], xu[i], yu[i], zu[i], q[i], aa, bb, cc = map(float, inp.readline().split())
           if atype[i] == 1: ncs += 1
           if atype[i] == 4: ncs += 1
           if atype[i] == 2: npb += 1
           if atype[i] == 3: ni += 1
           if atype[i] == 2: xpb.append(xu[i])
           if atype[i] == 2: indpb.append(i)
       qtot += q[i]
       ncube = int(np.power(npb+0.1,1./3.))
       if nframe == 1:
          print ('cube # = %i' %ncube)
          print ('total charge = %5.3lf' %qtot)
          print ('')
       if npb != len(indpb): print ('wrong Pb #', npb, len(indpb))
       line = inp.readline()

# -- Done reading trajectory --

# move system c.o.m. to the origin
xcom = np.mean(xu)
ycom = np.mean(yu)
zcom = np.mean(zu)
xu -= xcom
yu -= ycom
zu -= zcom
if abs(np.mean(xu))>0.001: print ('x_c.o.m. not zero?')
#print ('x c.o.m. = %lf  ' %(np.mean(xu)))

# -- define surface and lattice constant --
xmax = -10000.0
ymax = -10000.0
zmax = -10000.0
xmin = 10000.0
ymin = 10000.0
zmin = 10000.0
for j in range(int(ntot)):
    if xu[j] > xmax: xmax = xu[j]
    if yu[j] > ymax: ymax = yu[j]
    if zu[j] > zmax: zmax = zu[j]
    if xu[j] < xmin: xmin = xu[j]
    if yu[j] < ymin: ymin = yu[j]
    if zu[j] < zmin: zmin = zu[j]
ins = xmax-xmin
latpb = (max(xpb)-min(xpb))/float(ncube-1.)
latc = ins/float(ncube)
rcri = latpb*float(ncube)/2.0 - latpb*0.3
print ('system length = %lf' %(ins))
print ('lattice constant = %lf ' %latc)
print ('lattice constant from Pb = %lf ' %latpb)
print ('Surface at %lf and %lf' %(xmax, xmin))
print ('')

# -- check surface I atoms --
nis = 0 # of surface I atoms
rmv = []
for j in range(int(ntot)):
    if (abs(xu[j])>rcri or abs(yu[j])>rcri or abs(zu[j])>rcri) and atype[j] == 3: nis += 1
    if (abs(xu[j])>rcri or abs(yu[j])>rcri or abs(zu[j])>rcri) and atype[j] == 3: rmv.append(j)
nistheory = float(ncube*ncube)*6
if nis != nistheory: print ('wrong # of surface I', nis, nistheory)

indi = []
for j in range(int(ntot)):
    if j not in rmv and atype[j] == 3: indi.append(j)
if ni != len(indi) + nis: print ('wrong I #', ni, len(indi), nis)


# -- For each I, find two Pb index for each I, compute angle --
angle = []
walld = []
for i in range(len(indi)):
    nchk = []
    ii = indi[i] # index of iodide
    if atype[ii] != 3: print ('!!! Not I ? !!!')
    if i % 200 == 0: print ('%i th I atom out of %i' %(i, len(indi)))
    alist = [ abs(xmax-xu[ii]), abs(xu[ii]-xmin), abs(ymax-yu[ii]), abs(yu[ii]-ymin), abs(zmax-zu[ii]), abs(zu[ii]-zmin) ]
    walld.append(min(alist))

    # find two Pb index
    for j in range(npb):
        jpb = indpb[j]
        if atype[jpb] != 2: print ('!!! Not Pb ? !!!')
        dx = xu[ii] - xu[jpb]
        dy = yu[ii] - yu[jpb]
        dz = zu[ii] - zu[jpb]
        rr = np.power(dx*dx + dy*dy + dz*dz, 0.5) # PbI distance
        if rr < latpb*0.9: nchk.append(jpb)
    if len(nchk) != 2: print ('# of neighboring Pb is not 2?', len(nchk))

    # compute angle
    ipb1 = nchk[0]
    ipb2 = nchk[1]
    dx1 = xu[ipb1] - xu[ii]	# 1st Pb-I vector
    dy1 = yu[ipb1] - yu[ii]
    dz1 = zu[ipb1] - zu[ii]
    dx2 = xu[ipb2] - xu[ii]	# 2nd Pb-I vector
    dy2 = yu[ipb2] - yu[ii]
    dz2 = zu[ipb2] - zu[ii]
    rr1 = np.power(dx1*dx1 + dy1*dy1 + dz1*dz1, 0.5)
    rr2 = np.power(dx2*dx2 + dy2*dy2 + dz2*dz2, 0.5)
    inner = dx1*dx2 + dy1*dy2 + dz1*dz2
    acos = inner / rr1 / rr2
    if acos > 1.0: acos -= np.power(10.,-7.)
    if acos < -1.0: acos += np.power(10.,-7.)
    atheta = math.degrees(math.acos(acos))
    angle.append(atheta)

if len(walld) != len(indi): print ('something wrong in walld', len(walld), len(indi))


# -- Make histogram for avg. angle --
drr = 2.0
print ('')
print ('Avg. histogram bin size = %lf' %drr)
print ('Avg. PbIPb angle = %lf' %(np.mean(angle)))
ndr = int(max(angle)/drr) + int(5.0/drr)
hist = np.zeros(ndr)
for i in range(len(angle)):
    kk = int(angle[i]/drr) 
    hist[kk] += 1.0

# -- Make histogram for avg. angle vs. surf_dist --
dsr = 2.0
print ('Avg. surface hist. bin size = %lf' %dsr)
nsr = int(max(walld)/dsr) + int(3.0/dsr)
ncnt = np.zeros(nsr)
shist = np.zeros(nsr)
shist2 = np.zeros(nsr)
for i in range(len(angle)):
    ks = int(walld[i]/dsr) 
    ncnt[ks] += 1.0
    shist[ks] += angle[i]
    shist2[ks] += angle[i]*angle[i]

# -- Printings --
print ('')
output1.write('# PbIPb ang.   prob.   \n')
probtot = 0.0
for i in range(ndr):
    iang = float(i) * drr
    prob = hist[i]/float(len(angle))
    probtot += prob
    output1.write('%lf   %lf   \n' %(iang, prob))
print ('Avg. PbI dist prob. tot = %lf ' %probtot)

output2.write('# dist. from surf.    avg. PbIPb ang.  \n')
for j in range(nsr):
    jdist = float(j) * dsr
    if ncnt[j] > 0.0: 
       avgd = shist[j]/float(ncnt[j])
       avgd2 = shist2[j]/float(ncnt[j])
       varr = avgd2 - avgd*avgd
       if abs(varr) < np.power(10.,-10.): stdd = 0.0
       if abs(varr) >= np.power(10.,-10.): stdd = np.power(varr,0.5)
       output2.write('%lf   %lf   %lf   \n' %(jdist, avgd, stdd))


print ('')
#print ('************************ total simulation time *************************')
#print ('tstart = %lf' %(tstart/1000.), 'ps', ' & tend = %lf' %(tend/1000.), 'ps')
#print ('')


