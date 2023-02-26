
import numpy as np
import os
import sys

filename = sys.argv[1]
if not os.path.isfile(filename):
   print ('File name {} does not exist.'.format(filename))
   sys.exit()

output1=open('plmp.dat', 'w')
output2=open('pinit.xyz', 'w')

print ('')

nunitx = 10    # of units in x-axis
lcut = 35.
qcs = 0.86
qpb = 0.85
qi = -0.57
atom_masses = [132.91, 207.20, 126.90]
atom_charges = [qcs, qpb, qi]

print ('charge: Pb %lf, Cs %lf, I %lf '%(qpb, qcs, qi))
print ('')
with open(filename) as inp:
   line = inp.readline()
   while len(line) != 0:
       blk = inp.readline()
       natom, dummy = map(str, inp.readline().split())
       ntot = float(natom)
       blk = inp.readline()
       ntype, dummy1, dummy2 = map(str, inp.readline().split())
       natom_types = int(ntype)
       blk = inp.readline()
       # read boundaries
       alow, ahi, dummy1, dummy2 = map(str, inp.readline().split())
       xlo = float(alow)
       xhi = float(ahi)
       alow, ahi, dummy1, dummy2 = map(str, inp.readline().split())
       ylo = float(alow)
       yhi = float(ahi)
       alow, ahi, dummy1, dummy2 = map(str, inp.readline().split())
       zlo = float(alow)
       zhi = float(ahi)
       blk = inp.readline()
       dummy = inp.readline()
       blk = inp.readline()
       # read Masses
       ai, amass, dm1, dm2 = map(str, inp.readline().split())
       if atom_masses[int(ai)-1] != float(amass): print ('wrong mass Cs')
       ai, amass, dm1, dm2 = map(str, inp.readline().split())
       if atom_masses[int(ai)-1] != float(amass): print ('wrong mass Pb')
       ai, amass, dm1, dm2 = map(str, inp.readline().split())
       if atom_masses[int(ai)-1] != float(amass): print ('wrong mass I')
       blk = inp.readline()
       dummy = inp.readline()
       blk = inp.readline()
       # read Atoms
       atomid = []
       molid = []
       atype = []
       q = []
       xpos = []
       ypos = []
       zpos = []
       xpb = []
       npb = 0
       ncs = 0
       ni = 0
       qtot = 0.0
       for i in range(int(ntot)):
           a1, a2, a3, aq, ax, ay, az = map(float, inp.readline().split())
           atomid.append(int(a1))
           molid.append(int(a2))
           atype.append(int(a3))
           q.append(aq)
           xpos.append(ax)
           ypos.append(ay)
           zpos.append(az)
           qtot += aq
           if int(a3) == 1: ncs += 1
           if int(a3) == 3: ni += 1
           if int(a3) == 2: npb += 1
           if int(a3) == 2: xpb.append(ax)
       if len(xpos) != int(ntot): print ('wrong total # of atoms')
       ncube = int(np.power(npb+0.1,1./3.))
       print ('cube # = %i' %ncube)
       print ('total charge = %5.3lf' %qtot)
       print ('')
       line = inp.readline()

# -- Done reading initial for all Cs --

# move system c.o.m. to the origin
xcom = np.mean(xpos)
ycom = np.mean(ypos)
zcom = np.mean(zpos)
xpos -= xcom
ypos -= ycom
zpos -= zcom
#print ('x c.o.m. = %lf  ' %(np.mean(xpos)))

# -- define lattice constant --
xmax = -10000.0
xmin = 10000.0
for j in range(int(ntot)):
    if xpos[j] > xmax: xmax = xpos[j]
    if xpos[j] < xmin: xmin = xpos[j]
ins = xmax-xmin
latpb = (max(xpb)-min(xpb))/float(ncube-1.)
latc = ins/float(ncube)
rcri = latpb*float(ncube)/2.0 - latpb*0.6
print ('x_max = %lf and x_min = %lf ' %(xmax, xmin))
print ('system length = %lf' %(ins))
print ('lattice constant = %lf ' %latc)
print ('lattice constant from Pb = %lf ' %latpb)
print ('surface rcut = %lf ' %rcri)
print ('')

# -- define Corner Cs atoms --
ncnr = 0 # of corner Cs atoms
rmv = []
for j in range(int(ntot)):
    if (abs(xpos[j])>rcri and abs(ypos[j])>rcri and abs(zpos[j])>rcri) and atype[j] == 1: ncnr += 1
    if (abs(xpos[j])>rcri and abs(ypos[j])>rcri and abs(zpos[j])>rcri) and atype[j] == 1: rmv.append(j)
print ('# of corner Cs atoms = %i ' %ncnr)

# -- update all without corner Cs
newx = []
newy = []
newz = []   # new z-coodinate
newq = []   # new charge
newt = []   # new atom type
nqtot = 0.0
for j in range(int(ntot)):
    if j not in rmv:
       newx.append(xpos[j])
       newy.append(ypos[j])
       newz.append(zpos[j])
       newq.append(q[j])
       newt.append(atype[j])
       nqtot += q[j]
if len(newx)!= ntot-ncnr: print ('something wrong....!')
if len(newt)!= ntot-ncnr: print ('something wrong....!')
if len(newt)!= ntot-ncnr: print (len(newt), ntot, ncnr)
print ('total charge without corner Cs = %lf' %nqtot)
ichk = 0
for j in range(int(ntot)):
    if j in rmv: ichk += 1
    if j not in rmv:
       if xpos[j] != newx[j-ichk]: print ('x', ichk, j, xpos[j], newx[j-ichk])
       if ypos[j] != newy[j-ichk]: print ('y', ichk, j, ypos[j], newy[j-ichk])
       if zpos[j] != newz[j-ichk]: print ('z', ichk, j, zpos[j], newz[j-ichk])

# -- define surface Cs and assign new charge --
nsur = 0	# of surface Cs atoms
index = []
for j in range(int(ntot-ncnr)):
    if (abs(newx[j])>rcri or abs(newy[j])>rcri or abs(newz[j])>rcri) and newt[j] == 1: nsur += 1
    if (abs(newx[j])>rcri or abs(newy[j])>rcri or abs(newz[j])>rcri) and newt[j] == 1: index.append(j)
nstheory = (ncube-1)*(ncube-1)*6 + (ncube-1)*12
if nstheory != nsur: print ('check # of surface Cs again', nstheory)
print ('# of surface Cs atoms = %i ' %nsur)
qnow = float(npb)*qpb + float(ni)*qi + float(ncs-ncnr-nsur)*qcs
nqcs = -1.0*qnow/float(nsur)
print ('Surface Cs charge: %lf --> %lf ' %(qcs, nqcs))
print ('')

# -- assign new charge and type for surface Cs
for j in range(int(ntot-ncnr)):
    if j in index: newt[j] = 4
    if j in index: newq[j] = nqcs

# -- writing lammps new initial file
if natom_types+1!=max(newt): print (natom_types+1, max(newt))

atom_type_map = np.array(["Cs","Pb","I"])
atom_numbers = [55, 82, 53]
atom_masses = [132.91, 207.20, 126.90]
atom_charges = [qcs, qpb, qi]

output1.write('#LAMMPS configuration without high energy Cs \n')
output1.write('\n')
output1.write('%i atoms \n' %(ntot-ncnr))
output1.write('\n')
output1.write("%i atom types\n" %(natom_types+1))
output1.write("\n")
output1.write("%lf   %lf    xlo xhi\n" %(xlo, xhi))
output1.write("%lf   %lf    ylo yhi\n" %(ylo, yhi))
output1.write("%lf   %lf    zlo zhi\n" %(zlo, zhi))
output1.write("\n")
output1.write("Masses\n")
output1.write("\n")
# Masses
for i in range(natom_types):
    output1.write("%i   %lf   # %s\n" %(i+1, atom_masses[i], atom_type_map[i] ))
# Add additional Cs at the end
output1.write("%i   %lf   # %s\n" %(natom_types+1, atom_masses[0], atom_type_map[0] ))
output1.write("\n")
output1.write("Atoms\n")
output1.write("\n")
# Atoms
mqtot = 0.0
for i in range(int(ntot-ncnr)):
    mqtot += newq[i]
    output1.write('%i   %i   %i   %lf   %lf   %lf   %lf   \n' %(i+1, i+1, newt[i], newq[i], newx[i], newy[i], newz[i]))
print ('Total charge check --> %lf' %mqtot)

# writing vmd file
output2.write('%i \n' %(int(ntot-ncnr)))
output2.write('Atoms \n')
for k in range(int(ntot-ncnr)):
    if int(newt[k]) > natom_types+1: print ('check atom type')
    if int(newt[k]) == 1: output2.write('Cs ')
    if int(newt[k]) == 2: output2.write('Pb ')
    if int(newt[k]) in [3]: output2.write('I ')
    if int(newt[k]) in [4]: output2.write('Br ') # actually surface Cs
    output2.write('  %lf' %(newx[k] + 0.0))
    output2.write('  %lf' %(newy[k] + 0.0))
    output2.write('  %lf \n' %(newz[k] + 0.0))

print ('')
#print ('************************ total simulation time *************************')
#print ('tstart = %lf' %(tstart/1000.), 'ps', ' & tend = %lf' %(tend/1000.), 'ps')
#print ('')


