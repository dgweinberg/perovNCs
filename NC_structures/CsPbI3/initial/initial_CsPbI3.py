
import numpy as np
import os

## This initial file is for one layer

# Parameters ----------------------------------------------
# Pb --> 1 / I --> 2 / Cs --> 3

# Masses
m = np.zeros(3)
m[0] = 207.2		# Pb
m[1] = 126.90     # I 
m[2] = 132.91		# Cs

# charges 
q1 = np.zeros(3)
q1[0] = 0.85
q1[1] = -0.57
q1[2] = 0.86

# number of layers in each axis
nunity = 1     # number of unit in y axis
nunitx = 8     # number of inorganic part in x axis
nunitz = 8     # number of inorganic part in z axis
nlayer = 8     # number of inorganic layers in a unit cell

# length of unit cell
dist = 3.14470*1.0    # distance btwn Pb and I
yinorg = dist*2.0
xinorg = dist*2.0
zinorg = dist*2.0
yorg = 3.14470*1.0  #5.6 + 0.5

# length of whole repeated unit layer
ylen = dist*2.0 #yinorg*nunity + yorg*2.0 + 1.5  # length of y-unit layer
xlen = xinorg*nunitx     # length of x-unit layer
zlen = zinorg*nunitz     # length of z-unit layer

# number of atoms
nin = 4
no = 1                 # number of atoms in a ligand
nuin = 0 + nin*nunity   # number of inorg. atoms in a unit cell
nuo = no*1              # number of org. atoms in a unit cell
nuc = nuin + nuo        # number of atoms in a unit cell
ninorg = nin*nunitx*nunity*nunitz #+ nunitx*nunitz
norg = (no*nunitx*nunitz)*1
ntot = (ninorg + norg)*nlayer # total number of atoms

natype = 3 #no + 3         # number of atom types
numol = 5               # number of molecules in a unit cell

nblig = 0              # number of bonds in a ligand
nbuc = nblig*2          # number of bonds in a unit cell (w/ two ligands)
nbtype = 0              # number of bond types in a unit cell

nalig = 0               # number of angles in a ligand
nauc = nalig*2          # number of angles in a unit cell
ngtype = 0              # number of angle types

ndlig = 0               # number of dihedrals in a ligand
nduc = ndlig*2          # number of dihedrals in a unit cell
ndtype = 0              # number of dihedral types

nmol = numol*nunitx*nunitz*nlayer      # total number of molecules
nbtot = nbuc*nunitx*nunitz*nlayer      # total number of bonds
ngtot = nauc*nunitx*nunitz*nlayer      # total number of angles
ndtot = nduc*nunitx*nunitz*nlayer      # total number of dihedrals

# ---------------------------------------------------------
# Coordinates of a unit cell ------------------------------
# inorganic part
#inorg = np.loadtxt('PbX')
#ibelow = inorg[-1,:]
#iunit = inorg[:4,:]

iunit = np.zeros((4,4))
iunit[0,:] = [1, 0.0, 0.0, 0.0]		# Pb
iunit[1,:] = [2, 0.0, dist, 0.0]		# I
iunit[2,:] = [2, 0.0, 0.0, dist]		# I
iunit[3,:] = [2, dist, 0.0, 0.0]		# I
#ibelow = [1, 0.0, -dist, 0.0]		# I

#iyunit = ibelow
#inref = ibelow

iyadd = iunit
for i in range(nunity):
    iyadd[:,2] = iunit[:,2] + yinorg*i
    iyunit = iyadd #np.vstack( [iyunit, iyadd] )
    inref = iyadd #np.vstack( [inref, iyadd] )

# down ligand
#ligdown = np.loadtxt('down-nBA.txt')
#ldref = np.loadtxt('down-nBA.txt')

#ligdown[:,2] = ligdown[:,2] - dist - 0.5     # shift in y-axis
#ligdown[:,1] = ligdown[:,1] + dist           # shift in x-axis
#ligdown[:,3] = ligdown[:,3] + dist           # shift in z-axis
#ldref[:,2] = ldref[:,2] - dist - 0.5
#ldref[:,1] = ldref[:,1] + dist
#ldref[:,3] = ldref[:,3] + dist

# up ligand
ligup = np.loadtxt('up-Cs.txt')
luref = np.loadtxt('up-Cs.txt')

ligup[2] = ligup[2] + dist
ligup[1] = ligup[1] + dist
ligup[3] = ligup[3] + dist
luref[2] = luref[2] + dist
luref[1] = luref[1] + dist
luref[3] = luref[3] + dist

unit = iyunit
unit = np.vstack( [unit, luref])
#unit = np.vstack( [unit, ldref])

# check error
if len(unit) != 0 + nin*nunity + 1*no:
   print ('check # of atoms in a layer')
   print (len(unit), 0 + nin*nunity + 1*no)

if nunity == 1 and nunitx == 1 and nunitz == 1:
   np.savetxt ('unit.txt', unit, fmt='%8.5lf')

# ---------------------------------------------------------
# span unit cell in x-z plane -----------------------------
uadd = iyunit
uadd = np.vstack( [uadd, ligup])
#uadd = np.vstack( [uadd, ligdown])
uaddx = iyunit
uaddx = np.vstack( [uaddx, ligup])
#uaddx = np.vstack( [uaddx, ligdown])

ixadd = uadd
for i in range(nunitx-1):
    ixadd[:,1] += xinorg
    uaddx = np.vstack( [uaddx, ixadd] )
    unit = np.vstack( [unit, ixadd] )

izadd = uaddx
for i in range(nunitz-1):
    izadd[:,3] += zinorg
    unit = np.vstack( [unit, izadd] )

# ---------------------------------------------------------
# Constructing whole big matrix for coordinates -----------

natom = ntot    # total number of atoms
coord = np.empty(shape=[0, 4])
cadd = unit
cadd[:,2] -= ylen
cadd[:,1] += 0.0 #dist
cadd[:,3] += 0.0 #dist
for i in range(nlayer):
    cadd[:,2] += ylen
    if i%2 == 0:
       cadd[:,1] -= 0.0 #dist
       cadd[:,3] -= 0.0 #dist
    elif i%2 == 1:
       cadd[:,1] += 0.0 #dist
       cadd[:,3] += 0.0 #dist
    coord = np.vstack( [coord, cadd] )
# coord = atom type / x / y / z

if len(coord) != ntot:
   print ('check total # atoms in a unit layer')
   print (ntot, ninorg, norg)
   print (coord)

# ---------------------------------------------------------
# Writing VMD xyz file -----------------------------------

output2=open('init.xyz','w')

output2.write('%i \n' %ntot)
output2.write('Atoms \n')
for i in range(len(coord)):
    if coord[i,0] == 1:
       output2.write('Pb')
    elif coord[i,0] in [2]:
       output2.write('I')
    elif coord[i,0] == 3:
       output2.write('Cs')
    output2.write('%9.5lf' %coord[i,1])
    output2.write('%9.5lf' %coord[i,2])
    output2.write('%9.5lf \n' %coord[i,3])

# ---------------------------------------------------------
# Writing map.in file (phonon calc.) ----------------------

output=open('map.in', 'w')

output.write('%i %i %i %i \n' %(nunitx, nlayer, nunitz, nuc))
output.write('#l1 l2 l3 k tag \n')

tag = [ a + 1 for a in range(ntot) ]

print ('ntot =', ntot)
print ('# of atoms', nuc * nunitx * nlayer * nunitz) 

# stacking sequence : x --> z --> y
for i in range(ntot):
    # x unit cell label
    nl1 = i // nuc
    l1 = nl1 % nunitx
    # y unit cell label
    nl2 = i // (nuc * nunitx * nunitz)
    l2 = nl2 % nlayer
    # z unit cell label
    nl3 = i // (nuc * nunitx)
    l3 = nl3 % nunitz
    # atom type & ID
    k = i % nuc
    t = tag[i]
    output.write('%i %i %i %i %i \n' %(l1, l2, l3, k, t))


# ---------------------------------------------------------
# Writing LAMMPS initial file -----------------------------

output1=open('init.dat','w')

# Header --------------------------------------------------
output1.write('LAMMPS INITIAL FILE \n')
output1.write('\n')
output1.write('%i  atoms \n' %ntot)
#output1.write('%i  bonds \n' %nbtot)
#output1.write('%i  angles \n' %ngtot)
#output1.write('%i  dihedrals \n' %ndtot)

#output1.write(' \n')
output1.write(' %i  atom types \n' %natype )
#output1.write(' %i  bond types \n' %nbtype )
#output1.write(' %i  angle types \n' %ngtype )
#output1.write(' %i  dihedral types \n' %ndtype )

ybox = ylen*nlayer
output1.write('\n')
output1.write('0.0   %9.4lf   xlo xhi \n' %xlen)
output1.write('0.0   %9.4lf   ylo yhi \n' %ybox)
output1.write('0.0   %9.4lf   zlo zhi \n' %zlen)
#output1.write('0.0  0.0  0.0    xy xz yz \n')

output1.write(' \n')
output1.write('Masses \n')
output1.write(' \n')
for i in range(natype): 
    output1.write('%d  %lf \n' %(i+1, m[i]))

# Atoms ---------------------------------------------------
# atom ID
arange = range(ntot)
atomid = [ a + 1 for a in arange ]

# molecule ID
molid = []
nmuc = 0
for i in range(ntot):
    molid.append(0.0)
    if i % nuc < nuin:
       nmuc += 1
       molid[-1] = nmuc
    elif nuin <= i % nuc < nuin+no:
       if i % nuc == nuin: nmuc += 1
       molid[-1] = nmuc
    elif i % nuc >= nuin+no:
       if i % nuc == nuin+no: nmuc += 1
       molid[-1] = nmuc

# check error
if nmuc != nmol:
   print ('check # of molcules')
   print (nmuc, nmol)

# atom type
atype = coord[:,0]

# charge
q = []   # charge array
for i in range(ntot):
    q.append(0.0)
    nn = int(coord[i,0])
    q[-1] = q1[nn-1]

# coordinates
x = coord[:,1]
y = coord[:,2]
z = coord[:,3]

# Writing
output1.write('\n')
output1.write('Atoms \n')
output1.write('\n')

# atom ID / molecule ID / atom type / charge / x-coord. / y-coord. / z-coord.
atoms = np.column_stack((atomid, molid, atype, q, x, y, z))
np.savetxt(output1, atoms, fmt='%d  %d  %d  %11.6lf  %11.6lf  %11.6lf  %11.6lf')

