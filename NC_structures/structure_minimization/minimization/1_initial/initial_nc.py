import numpy as np;
import math
import sys

#########################################################################################
#read inputs

if len(sys.argv)!=4 and len(sys.argv)!=5:
    print("Invalid number of arguments\nUsage: python perovMixBuild.py nx ny nz [par]\n")
    exit()


nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])


#par = 0 means cubic, par = 1 means ortho
if len(sys.argv)==5: 
	par = float(sys.argv[4])
	if(0>par or par>1):
		print("Need  0 <= par <= 1!")
		exit()
else: 
	par = 0.0

print ('')
if par == 0.0: print ('Cubic initial config.')
if par == 1.0: print ('Orthorhombic initial config.')
print ('')

if(nx*ny*nz>10000):
	print("Woah there! This dot is too big (nx: {nx} ny: {ny} nz: {nz})")
	exit()

#########################################################################################
#set up needed info
#########################################################################################
#set up cubic and ortho parameters


cubicScale = 11.884*math.sqrt(2)
cubicA = 1.0
cubicB = 1.0
cubicC = math.sqrt(2)

orthoScale = 16.736
orthoA = 1.0
orthoB = 0.96838
orthoC = 1.40831

#########################################################################################
#interpolate between

parScale = cubicScale*(1-par)+orthoScale*par
parA = cubicA*(1-par)+orthoA*par
parB = cubicB*(1-par)+orthoB*par
parC = cubicC*(1-par)+orthoC*par



d1 = (0.25-0.19462)*par
d2 = (0.25-0.19731)*par
d3 = 0.03577*par
d4 = 0.00113*par
d5 = 0.06202*par
d6 = 0.04005*par
d7 = 0.00509*par

#########################################################################################

#details for basic conf
#units in BOHR!!!

basis_vectors= parScale * np.array([[ parA/math.sqrt(2), -parA/math.sqrt(2), 0 ],[ parB/math.sqrt(2), parB/math.sqrt(2), 0 ], [ 0, 0, parC ]])
atoms = np.array(["Cs", "Cs", "Cs", "Cs", "Pb", "Pb", "Pb", "Pb", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I", "I"])

lxy = 0.5 * parScale*math.sqrt(parA**2+parB**2)
lz = 0.5 * parScale*parC

atom_positions = np.array([
[0.50-d6, 0.5+d7, 0.25000],
[0.50+d6, 0.5-d7, 0.75000],
[1.00-d6, 0.0+d7, 0.25000],
[0.00+d6, 1.0-d7, 0.75000],
[0.50000, 0.00000, 0.0000],
[0.00000, 0.50000, 0.0000],
[0.50000, 0.00000, 0.5000],
[0.00000, 0.50000, 0.5000],
[0.25-d1, 0.25-d2, 0.0+d3],
[0.75-d1, 0.25+d2, 0.0+d3],
[0.75+d1, 0.75+d2, 1.0-d3],
[0.25+d1, 0.75-d2, 1.0-d3],
[0.25-d1, 0.25-d2, 0.5-d3],
[0.75-d1, 0.25+d2, 0.5-d3],
[0.75+d1, 0.75+d2, 0.5+d3],
[0.25+d1, 0.75-d2, 0.5+d3],
[0.00+d4, 0.50+d5, 0.2500],
[0.50+d4, 1.00-d5, 0.2500],
[1.00-d4, 0.50-d5, 0.7500],
[0.50-d4, 0.00+d5, 0.7500]])


#details for lammps config
natom_types = 3
atom_type_map = np.array(["Cs","Pb","I"])
atom_numbers = [55, 82, 53]
atom_masses = [132.91, 207.20, 126.90]
atom_charges = [0.86, 0.85,-0.57]
formal_charges = [1,2,-1]


#########################################################################################
#build bulk to cut dot from
#cut the nanocrystal from bulk
xmin = -lxy * (math.floor(nx/2.0)+0.1)
xmax =  lxy * (math.ceil(nx/2.0)+0.1)

ymin = -lxy * (math.floor(ny/2.0)+0.1)
ymax =  lxy * (math.ceil(ny/2.0)+0.1)

zmin = -0.6*lz
zmax = lz * (nz-0.4)

#zmin = -0.5*lz
#zmax = lz * (nz-0.5)


built_dot = []


for x_cell in range(-nx,int(nx)+1):
	for y_cell in range (-ny, int(ny)+1):
		for z_cell in range(-1,int(nz)+1):
				for i,atom in enumerate(atoms):
					#pos = basis_vectors.transpose() @ atom_positions[i] + basis_vectors.transpose() @ [x_cell,y_cell,z_cell]
					pos = np.matmul( basis_vectors.transpose() , atom_positions[i] ) + np.matmul( basis_vectors.transpose() , [x_cell,y_cell,z_cell] )
					if pos[0]>=xmin and pos[0]<=xmax and pos[1]>=ymin and pos[1]<=ymax and pos[2] >= zmin and pos[2]<=zmax:
						#print(atom, pos[0], pos[1], pos[2])
						built_dot.append([atom, pos[0], pos[1], pos[2]])


#########################################################################################
#ballance formal charge

natoms = len(built_dot)

formal_charge = 0
surface_list = []
to_remove = []

#print(ny*basis_vectors[1][1])

for i in range(natoms):
	if built_dot[i][1]<= 0 or built_dot[i][2]<=0 or built_dot[i][3]<=0:
		surface_list.append(i)
	elif built_dot[i][1]>=(nx-1)*basis_vectors[0][0] or built_dot[i][2]>=(ny-1)*basis_vectors[1][1] or built_dot[i][3]>=(nz-1)*basis_vectors[2][2]:
		surface_list.append(i)

np.savetxt('surface.txt',surface_list, '%d')


#########################################################################################
#write outputs


#filter output


filter_of = open("conf.par", "w")
#filter_of.write(f"{natoms}\n")
filter_of.write("%i \n" %natoms)

for i in range(natoms):
    if not(i in to_remove):
      #filter_of.write(f"{built_dot[i][0]} {built_dot[i][1]} {built_dot[i][2]} {built_dot[i][3]}\n")
      filter_of.write("%s  %lf   %lf   %lf   \n" %(built_dot[i][0], built_dot[i][1], built_dot[i][2], built_dot[i][3]))
filter_of.close()


#lammps output
bohr_to_ang = 0.529177
built_dot = np.array(built_dot)
xpos=built_dot[:,1].astype(float)*bohr_to_ang
ypos=built_dot[:,2].astype(float)*bohr_to_ang
zpos=built_dot[:,3].astype(float)*bohr_to_ang

#print ('')
#print ('max x dist. in bohr = %lf' %((max(xpos)-min(xpos)) / bohr_to_ang))
#print ('max y dist. in bohr = %lf' %((max(ypos)-min(ypos)) / bohr_to_ang))
#print ('max z dist. in bohr = %lf' %((max(zpos)-min(zpos)) / bohr_to_ang))
print ('max x dist. in A = %lf' %((max(xpos)-min(xpos)) / 1.0))
print ('max y dist. in A = %lf' %((max(ypos)-min(ypos)) / 1.0))
print ('max z dist. in A = %lf' %((max(zpos)-min(zpos)) / 1.0))
print ('')

xcom = np.mean(xpos)
ycom = np.mean(ypos)
zcom = np.mean(zpos)
xpos -= xcom
ypos -= ycom
zpos -= zcom

xlo = min(xpos)-500
xhi = max(xpos)+500

ylo = min(ypos)-500
yhi = max(ypos)+500

zlo = min(zpos)-500
zhi = max(zpos)+500

total_charge = 0

lammps_of = open("lammpsconf.par", "w")
lammps_of.write("#LAMMPS configuration file for perovskite NanoCube (nx: {nx} ny: {ny} nz: {nz})\n")
lammps_of.write("\n")
#lammps_of.write(f"{natoms-formal_charge} atoms\n")
lammps_of.write("%i atoms\n" %(natoms-formal_charge))
lammps_of.write("\n")
lammps_of.write("%i atom types\n" %natom_types)
lammps_of.write("\n")
lammps_of.write("%lf   %lf    xlo xhi\n" %(xlo, xhi))
lammps_of.write("%lf   %lf    ylo yhi\n" %(ylo, yhi))
lammps_of.write("%lf   %lf    zlo zhi\n" %(zlo, zhi))
lammps_of.write("\n")
lammps_of.write("Masses\n")
lammps_of.write("\n")
for i in range(natom_types):
   lammps_of.write("%i   %lf   # %s\n" %(i+1, atom_masses[i], atom_type_map[i] ))
lammps_of.write("\n")
lammps_of.write("Atoms\n")
lammps_of.write("\n")

ncs = 0 # number of Cs atom
npb = 0
ni = 0
for i in range(natoms):
   atom_type = np.where(atom_type_map == built_dot[i,0])[0][0]
   if not(i in to_remove):
      total_charge+=atom_charges[atom_type]
      if atom_type+1 == 1: ncs += 1
      if atom_type+1 == 2: npb += 1
      if atom_type+1 == 3: ni += 1
      lammps_of.write("%i   %i   %i   %lf   %lf   %lf   %lf   \n" %(i+1, i+1, atom_type+1, atom_charges[atom_type], xpos[i], ypos[i], zpos[i]))
lammps_of.close()

#print ('')
print ('NC size = %lf' %(max(xpos)-min(xpos)))
print ('# of Cs = %i' %ncs)
print ('# of Pb = %i' %npb)
print ('# of I = %i' %ni)
qtot = 2.0*npb + 1.0*ncs - 1.0*ni
print("total charge: %lf " %total_charge)
print("formal charge: %lf " %qtot)
#nbnd = npb * 2 + ncs * 1 + ni * 17 # of valence electrons
#print("min # of states: %lf " %nbnd)
print ('')

#xyz output
xyz_of = open("initconf.xyz", "w")
xyz_of.write("%i \n" %natoms)
xyz_of.write("Atoms\n")
#xyz_of.write("#xyz file for perovskite NanoCube (nx: {nx} ny: {ny} nz: {nz})\n")
for i in range(natoms):
   if not(i in to_remove):
      xyz_of.write("%s   %lf   %lf   %lf   \n" %(built_dot[i,0], xpos[i], ypos[i], zpos[i]))
xyz_of.close()


exit()


