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
matplotlib.rcParams.update({'figure.figsize' : (8, 6),
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
                                    8*8*8,
                                    9*9*9,
                                    10*10*10, 11*11*11])**(1/3))


my_rect_size = Pb_to_nm(np.array([2*2*3, 2*3*3, 
                                     3*3*4, 3*4*4,
                                     4*4*5, 4*5*5, 
                                     5*5*6, 5*6*6, 
                                     6*6*7, 6*7*7,
                                     7*7*8, 7*8*8])**(1/3))


relxcubeFileDir = "../electronic/relax/"
relxrectFileDir = "../electronic/relax/"

cubiccubeFileDir = "../electronic/cubic/"
cubicrectFileDir = "../electronic/cubic/"

orthocubeFileDir = "../electronic/ortho/"
orthorectFileDir = "../electronic/ortho/"


fileExtension = "/bse/OS.dat"

cubeDirNames = ["2x2x2", "3x3x3", "4x4x4", "5x5x5", "6x6x6", "7x7x7", "8x8x8", "9x9x9", "10x10x10", "11x11x11"]
rectDirNames = ["2x2x3", "2x3x3", "3x3x4",
                "3x4x4", "4x4x5", "4x5x5", "5x5x6", "5x6x6",
                "6x6x7", "6x7x7", "7x7x8", "7x8x8"]
xsize = 5;
ysize = xsize*3/4
fig, ax = plt.subplots(figsize = (xsize,ysize))
sizes = []
DeltaE = []

for n,file in enumerate(cubeDirNames):
    print(f"Loading file: {relxcubeFileDir+file+fileExtension}")
    try:
        OS = np.loadtxt(relxcubeFileDir+file+fileExtension)
    except:
        print(f"Can't find: {relxcubeFileDir+file+fileExtension}")
    else:   
        if OS.size==0: continue; 
        OS[:,2]*=27211.396641308
        sizes.append(my_cubic_size[n])
        DeltaE.append((OS[1,2]+OS[2,2]+OS[3,2])/3.0-OS[0,2])

ax.plot(sizes,DeltaE, 'ko')

sizes = []
DeltaE = []

for n,file in enumerate(cubeDirNames):
    print(f"Loading file: {cubiccubeFileDir+file+fileExtension}")
    try:
        OS = np.loadtxt(cubiccubeFileDir+file+fileExtension)
    except:
        print(f"Can't find: {cubiccubeFileDir+file+fileExtension}")
    else:    
        if OS.size==0: continue;
        OS[:,2]*=27211.396641308
        sizes.append(my_cubic_size[n])
        DeltaE.append((OS[1,2]+OS[2,2]+OS[3,2])/3.0-OS[0,2])

ax.plot(sizes,DeltaE, 'bs')

sizes = []
DeltaE = []

for n,file in enumerate(cubeDirNames):
    print(f"Loading file: {orthocubeFileDir+file+fileExtension}")
    try:
        OS = np.loadtxt(orthocubeFileDir+file+fileExtension)
    except:
        print(f"Can't find: {orthocubeFileDir+file+fileExtension}")
    else: 
        if OS.size==0: continue;   
        OS[:,2]*=27211.396641308
        sizes.append(my_cubic_size[n])
        DeltaE.append((OS[1,2]+OS[2,2]+OS[3,2])/3.0-OS[0,2])

ax.plot(sizes,DeltaE, 'gd')



ax.set_ylabel("$\Delta$(meV)")
ax.set_xlabel("Size (nm)")
ax.set_ylim(ymax = 150, ymin = 0)
plt.xlim(xmax = 8, xmin = 1)
ax.set_title("Dark-Bright Splitting")

fig.savefig("Bright-Dark.png", dpi=300, bbox_inches='tight')