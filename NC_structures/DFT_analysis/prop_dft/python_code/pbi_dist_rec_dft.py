
import numpy as np
import os
import sys

filename = sys.argv[1]
index1 = sys.argv[2]
index2 = sys.argv[3]
index3 = sys.argv[4]
if not os.path.isfile(filename):
   print ('File name {} does not exist.'.format(filename))
   sys.exit()

# how to use example
# python pbi_dist_rec_dft.py prop.xyz 2 2 2

output1=open('avg_pbi.dat', 'w') # histogram of PbI dist.
output2=open('sur_pbi.dat', 'w') # surf_dist. vs. PbI dist.

nux = int(index1)
nuy = int(index2)
nuz = int(index3)
print ('')
print ('# of unitx in x-axis = %i' %nux)
print ('# of unitx in y-axis = %i' %nuy)
print ('# of unitx in z-axis = %i' %nuz)
print ('')

dt = 1.0
nframe = 0     # of frames (snapshots)
#print ('')
#print ('# compute distribution for PbI distance #')
#print ('')

with open(filename) as inp:
   line = 'ing' #inp.readline()
   while len(line) != 0:
       # converting lammps dump file to VMD readable xyz format
       nframe += 1
       nch = 0
       #tstep = inp.readline()
       #line = inp.readline()
       #if nframe % 100 == 1: print ('%i th frame' %nframe)
       #if nframe == 1:
       #   tstart = int(tstep)*dt
       #tend = int(tstep)*dt
       natom = inp.readline()
       line = inp.readline()
       #xlo, xhi = map(float, inp.readline().split())
       #ylo, yhi = map(float, inp.readline().split())
       #zlo, zhi = map(float, inp.readline().split())
       # ----------------------------
       ntot = int(natom)
       #line = inp.readline()
       atype = np.zeros(ntot)
       xu = np.zeros(ntot)
       yu = np.zeros(ntot)
       zu = np.zeros(ntot)
       #q = np.zeros(ntot)
       nout = 0   # of outside atoms
       npb = 0    # of pb atoms
       ncs = 0
       ni = 0
       qtot = 0.0
       xpb = []
       ypb = []
       zpb = []
       indpb = []
       indi = []
       for i in range(ntot):
           atype[i], xu[i], yu[i], zu[i] = map(float, inp.readline().split())
           if atype[i] == 1: ncs += 1
           if atype[i] == 4: ncs += 1
           if atype[i] == 2: npb += 1
           if atype[i] == 3: ni += 1
           if atype[i] == 2: xpb.append(xu[i])
           if atype[i] == 2: ypb.append(yu[i])
           if atype[i] == 2: zpb.append(zu[i])
           if atype[i] == 2: indpb.append(i)
           if atype[i] == 3: indi.append(i)
       qtot += 0.0 #q[i]
       #ncube = int(np.power(npb+0.1,1./3.))
       #if nframe == 1:
          #print ('cube # = %i' %ncube)
          #print ('total charge = %5.3lf' %qtot)
          #print ('')
       if npb != len(indpb): print ('wrong Pb #', npb, len(indpb))
       if ni != len(indi): print ('wrong I #', ni, len(indi))
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
latpbx = (max(xpb)-min(xpb))/float(nux-1.)
latpby = (max(ypb)-min(ypb))/float(nuy-1.)
latpbz = (max(zpb)-min(zpb))/float(nuz-1.)
latpb = min(latpbx, latpby, latpbz)
latcx = ins/float(nux)
rcrix = latpbx*float(nux)/2.0 - latpbx*0.2
rcriy = latpby*float(nuy)/2.0 - latpby*0.2
rcriz = latpbz*float(nuz)/2.0 - latpbz*0.2
print ('system length = %lf' %(ins))
print ('lattice constant = %lf ' %latcx)
print ('lattice constant from Pb (x) = %lf ' %latpbx)
print ('Surface at %lf and %lf' %(xmax, xmin))
print ('Surface criteria at %lf' %(rcrix))
print ('')

# -- check Corner Cs atoms (shouldn't be here) --
ncnr = 0 # of corner Cs atoms
rmv = []
for j in range(int(ntot)):
    if (abs(xu[j])>rcrix and abs(yu[j])>rcriy and abs(zu[j])>rcriz) and atype[j] == 1: ncnr += 1
    if (abs(xu[j])>rcrix and abs(yu[j])>rcriy and abs(zu[j])>rcriz) and atype[j] == 1: rmv.append(j)
#if abs(ncnr)>0.01: print ('there is corner Cs atom..?')
print ('# of corner Cs atoms = %i ' %ncnr)


# -- For each Pb-I pair, check distance --
pbid = []
walld = []
for i in range(npb):
    npair = 0.0
    ipb = indpb[i]
    if atype[ipb] != 2: print ('!!! Not Pb ? !!!')
    #if i % 100 == 0: print ('%i th Pb atom' %i)
    for j in range(ni):
        ji = indi[j]
        if atype[ji] != 3: print ('!!! Not iodide ? !!!')
        dx = xu[ipb] - xu[ji]
        dy = yu[ipb] - yu[ji]
        dz = zu[ipb] - zu[ji]
        rr = np.power(dx*dx + dy*dy + dz*dz, 0.5) # PbI distance
        alist = [ abs(xmax-xu[ipb]), abs(xu[ipb]-xmin), abs(ymax-yu[ipb]), abs(yu[ipb]-ymin), abs(zmax-zu[ipb]), abs(zu[ipb]-zmin) ]
        if rr < latpb*1.0: pbid.append(rr)
        if rr < latpb*1.0: npair += 1.0
        if rr < latpb*1.0: walld.append(min(alist)) # min. dist. btwn Pb and surf.
    if npair != 6: print ('# of neighboring I is not 6?', npair)
if len(pbid) != npb*6: print ('check # of PbI pair', len(pbid), npb)
if len(walld) != len(pbid): print ('something wrong in walld', len(walld), len(pbid))

print ('# of Pb atoms = %i' %npb)
print ('# of I atoms = %i' %ni)
print ('# of Cs atoms = %i' %ncs)
print ('# of bonds = %lf' %(len(pbid)))

# -- Make histogram for avg. PbI dist. --
drr = 0.02
print ('')
print ('Avg. histogram bin size = %lf' %drr)
print ('Avg. PbI distance = %lf' %(np.mean(pbid)))
ndr = int(max(pbid)/drr) + int(1.0/drr)
hist = np.zeros(ndr)
for i in range(len(pbid)):
    kk = int(pbid[i]/drr) 
    hist[kk] += 1.0

# -- Make histogram for avg. PbI dist. vs. surf_dist --
dsr = 0.5
print ('Avg. surface hist. bin size = %lf' %dsr)
nsr = int(max(walld)/dsr) + int(3.0/dsr)
ncnt = np.zeros(nsr)
shist = np.zeros(nsr)
shist2 = np.zeros(nsr)
for i in range(len(pbid)):
    ks = int(walld[i]/dsr) 
    ncnt[ks] += 1.0
    shist[ks] += pbid[i]
    shist2[ks] += pbid[i]*pbid[i]

# -- Printings --
print ('')
output1.write('# PbI dist.   prob.   \n')
probtot = 0.0
for i in range(ndr):
    idist = float(i) * drr
    prob = hist[i]/float(len(pbid))
    probtot += prob
    output1.write('%lf   %lf   \n' %(idist, prob))
print ('Avg. PbI dist prob. tot = %lf ' %probtot)

output2.write('# dist. from surf.    avg. PbI dist.  \n')
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


