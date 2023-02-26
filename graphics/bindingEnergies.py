import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def nm_to_Pb(nm):
    Pb = nm/0.628874;
    return Pb
def Pb_to_nm(Pb):
    nm = Pb*0.628874;
    return nm;

# tick and axes widths, etc
matplotlib.rcParams.update({'figure.figsize' : (4, 3),
    'axes.linewidth' : 2.0,
    'xtick.major.size' : 3.0,
    'xtick.major.width' : 1.5,
    'ytick.major.size' : 3.0,
    'ytick.major.width' : 1.5,
    'font.size' : 16,
    'legend.fontsize' : 14,
    'xtick.labelsize' : 14,
    'ytick.labelsize' : 14,
    'text.usetex': True,
    'font.family': 'sans-serif',
    #'font.serif': 'Palotino',
    'font.weight': 'bold'})
my_cubic_size = Pb_to_nm(np.array([2*2*2, 
                                    3*3*3,
                                    4*4*4,
                                    5*5*5,
                                    6*6*6,
                                    7*7*7,
                                    8*8*8, 9*9*9, 10*10*10, 11*11*11])**(1/3))


my_rect_size = Pb_to_nm(np.array([2*2*3, 2*3*3, 
                                     3*3*4, 3*4*4,
                                     4*4*5, 4*5*5, 
                                     5*5*6, 5*6*6, 
                                     6*6*7, 6*7*7,
                                     7*7*8, 7*8*8])**(1/3))


cubeFileDir = "../electronic/cubic/"
rectFileDir = "../electronic/cubic/"
fileExtension = "/bse/exciton.dat"

cubeDirNames = ["2x2x2", "3x3x3", "4x4x4", "5x5x5", "6x6x6", "7x7x7", "8x8x8", "9x9x9"]#, "10x10x10", "11x11x11"]
rectDirNames = ["2x2x3", "2x3x3", "3x3x4",
                "3x4x4", "4x4x5", "4x5x5", "5x5x6", "5x6x6",
                "6x6x7", "6x7x7", "7x7x8", "7x8x8"]

fig, ax = plt.subplots(figsize = (4,3))

for n,file in enumerate(cubeDirNames):
    print(f"Loading file: {cubeFileDir+file+fileExtension}")
    try:
        OS = np.loadtxt(cubeFileDir+file+fileExtension)[1:4,6]
    except:
        print(f"Can't find: {cubeFileDir+file+fileExtension}")
    else:
        if (np.size(OS)==0): continue    
        ax.scatter(my_cubic_size[n]*np.ones_like(OS[0]), -1*OS[0], c = "Black", marker = 's', alpha=1)
        ax.scatter(my_cubic_size[n]*np.ones_like(OS[1:3]), -1*OS[1:3], c ='Red', marker = 's', alpha=1)


#plt.legend()
plt.xlim(xmax = 10, xmin = 1)
plt.ylabel("Exciton Binding Energy (eV)")
plt.xlabel("Size (nm)")
plt.xlim(xmax = 10, xmin = 1)
plt.title("Cubic Structures")
fig.savefig('CubicBindingEnergies.png', dpi=300, bbox_inches='tight')

###########################################

cubeFileDir = "../electronic/ortho/"
rectFileDir = "../electronic/ortho/"

fig, ax = plt.subplots(figsize = (4,3))

for n,file in enumerate(cubeDirNames):
    print(f"Loading file: {cubeFileDir+file+fileExtension}")
    try:
        OS = np.loadtxt(cubeFileDir+file+fileExtension)[1:4,6]
    except:
        print(f"Can't find: {cubeFileDir+file+fileExtension}")
    else:
        if (np.size(OS)==0): continue    
        ax.scatter(my_cubic_size[n]*np.ones_like(OS[0]), -1*OS[0], c = "Black", marker = 's', alpha=1)
        ax.scatter(my_cubic_size[n]*np.ones_like(OS[1:3]), -1*OS[1:3], c ='Red', marker = 's', alpha=1)


#plt.legend()
plt.xlim(xmax = 10, xmin = 1)
plt.ylabel("Exciton Binding Energy (eV)")
plt.xlabel("Size (nm)")
plt.xlim(xmax = 10, xmin = 1)
plt.title("Orthorhombic Structures")
fig.savefig('OrthoBindingEnergies.png', dpi=300, bbox_inches='tight')

###########################################

cubeFileDir = "../electronic/relax/"
rectFileDir = "../electronic/relax/"

fig, ax = plt.subplots(figsize = (4,3))

for n,file in enumerate(cubeDirNames):
    print(f"Loading file: {cubeFileDir+file+fileExtension}")
    try:
        OS = np.loadtxt(cubeFileDir+file+fileExtension)[1:4,6]
    except:
        print(f"Can't find: {cubeFileDir+file+fileExtension}")
    else:
        if (np.size(OS)==0): continue    
        ax.scatter(my_cubic_size[n]*np.ones_like(OS[0]), -1*OS[0], c = "Black", marker = 'o', alpha=1)
        #ax.scatter(my_cubic_size[n]*np.ones_like(OS[1:3]), -1*OS[1:3], c ='Red', marker = 'o', alpha=1)


#plt.legend()
plt.xlim(xmax = 10, xmin = 1)
plt.ylabel("Energy (eV)")
plt.xlabel("Size (nm)")
plt.xlim(xmax = 10, xmin = 1)
plt.title("Exciton Binding Energy")
fig.savefig('RelaxBindingEnergies.png', dpi=300, bbox_inches='tight')