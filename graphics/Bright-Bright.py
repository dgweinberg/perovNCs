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


prop_cycle = plt.rcParams['axes.prop_cycle']
# colors = prop_cycle.by_key()['color']
colors = ["grey", "grey", "grey", "grey", "grey", "grey", "grey"]

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


# https://doi.org/10.1038/s41563-022-01349-4
han_sizes = np.array([4.9,5.4,6.2,6.7,7.4,7.9])
han_errorbars = np.array([0.2,0.3,0.3,0.45,0.6,0.8])
han_TT_Splittings = np.array([1.61,1.4,1.25,1.20,0.85,0.7])


# https://link.aps.org/doi/10.1103/PhysRevLett.119.026401
yin_sizes = np.array([9.3])
yin_splittings = np.array([0.40])


# https://link.aps.org/doi/10.1103/PhysRevB.97.235304
nestoklon_sizes = np.array([10.0])
nestoklon_splittings = np.array([0.12])


#https://doi.org/10.1063/5.0027972
rossi2020_sizes=np.array([4,6.5])
rossi2020_splittings = np.array([7.5,2])


relxcubeFileDir = "/home/dweinberg/Projects/rashba/calcs/production2/relax/"
relxrectFileDir = "/home/dweinberg/Projects/rashba/calcs/production2/relax/"

cubiccubeFileDir = "/home/dweinberg/Projects/rashba/calcs/production2/cubic/"
cubicrectFileDir = "/home/dweinberg/Projects/rashba/calcs/production2/cubic/"

orthocubeFileDir = "/home/dweinberg/Projects/rashba/calcs/production2/ortho/"
orthorectFileDir = "/home/dweinberg/Projects/rashba/calcs/production2/ortho/"


fileExtension = "/bse/OS.dat"

cubeDirNames = ["2x2x2", "3x3x3", "4x4x4", "5x5x5", "6x6x6", "7x7x7", "8x8x8", "9x9x9"]#, "10x10x10", "11x11x11"]
rectDirNames = ["2x2x3", "2x3x3", "3x3x4",
                "3x4x4", "4x4x5", "4x5x5", "5x5x6", "5x6x6",
                "6x6x7", "6x7x7", "7x7x8", "7x8x8"]


xsize = 4;
ysize = xsize*3/4
fig, (axOrtho, axRelax) = plt.subplots(2,1,figsize = (xsize,ysize),sharex=True)

axRelax.errorbar(han_sizes,  han_TT_Splittings/2.0,xerr=han_errorbars, fmt = '.',mfc=colors[0], mec=colors[0], ecolor=colors[0], zorder=4)
axRelax.errorbar(han_sizes, -han_TT_Splittings/2.0,xerr=han_errorbars, fmt = '.',mfc=colors[0], mec=colors[0], ecolor=colors[0], zorder=4)
axOrtho.errorbar(han_sizes,  han_TT_Splittings/2.0,xerr=han_errorbars, fmt = '.',mfc=colors[0], mec=colors[0], ecolor=colors[0], zorder=4)
axOrtho.errorbar(han_sizes, -han_TT_Splittings/2.0,xerr=han_errorbars, fmt = '.',mfc=colors[0], mec=colors[0], ecolor=colors[0], zorder=4)


axRelax.errorbar(rossi2020_sizes,  rossi2020_splittings/2.0, fmt = '.',mfc=colors[1], mec=colors[1], ecolor=colors[1], zorder=3)
axRelax.errorbar(rossi2020_sizes, -rossi2020_splittings/2.0, fmt = '.',mfc=colors[1], mec=colors[1], ecolor=colors[1], zorder=3)
axOrtho.errorbar(rossi2020_sizes,  rossi2020_splittings/2.0, fmt = '.',mfc=colors[1], mec=colors[1], ecolor=colors[1], zorder=3)
axOrtho.errorbar(rossi2020_sizes, -rossi2020_splittings/2.0, fmt = '.',mfc=colors[1], mec=colors[1], ecolor=colors[1], zorder=3)


axRelax.errorbar(yin_sizes,  yin_splittings/2.0, fmt = '.',mfc=colors[2], mec=colors[2], ecolor=colors[2], zorder=2)
axRelax.errorbar(yin_sizes, -yin_splittings/2.0, fmt = '.',mfc=colors[2], mec=colors[2], ecolor=colors[2], zorder=2)
axOrtho.errorbar(yin_sizes,  yin_splittings/2.0, fmt = '.',mfc=colors[2], mec=colors[2], ecolor=colors[2], zorder=2)
axOrtho.errorbar(yin_sizes, -yin_splittings/2.0, fmt = '.',mfc=colors[2], mec=colors[2], ecolor=colors[2], zorder=2)

axRelax.errorbar(nestoklon_sizes, -nestoklon_splittings/2.0, fmt = '.',mfc=colors[3], mec=colors[3], ecolor=colors[3], zorder=1)
axRelax.errorbar(nestoklon_sizes,  nestoklon_splittings/2.0, fmt = '.',mfc=colors[3], mec=colors[3], ecolor=colors[3], zorder=1)
axOrtho.errorbar(nestoklon_sizes, -nestoklon_splittings/2.0, fmt = '.',mfc=colors[3], mec=colors[3], ecolor=colors[3], zorder=1)
axOrtho.errorbar(nestoklon_sizes,  nestoklon_splittings/2.0, fmt = '.',mfc=colors[3], mec=colors[3], ecolor=colors[3], zorder=1)

sizes = []
DeltaE1 = []
DeltaE2 = []


sizes = []
DeltaE1 = []
DeltaE2 = []

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
        DeltaE1.append(OS[2,2]-OS[1,2])
        DeltaE2.append(OS[3,2]-OS[1,2])
DeltaE1 = np.array(DeltaE1)
DeltaE2 = np.array(DeltaE2)
ave = (DeltaE1+DeltaE2)/3.0
axOrtho.plot(sizes,-ave, 'gd', zorder=10)
axOrtho.plot(sizes,DeltaE1-ave, 'gd', zorder=10)
axOrtho.plot(sizes,DeltaE2-ave, 'gd', zorder=10)
axOrtho.set_ylim(-4,4)


sizes = []
DeltaE1 = []
DeltaE2 = []

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
        DeltaE1.append(OS[2,2]-OS[1,2])
        DeltaE2.append(OS[3,2]-OS[1,2])
DeltaE1 = np.array(DeltaE1)
DeltaE2 = np.array(DeltaE2)
ave = (DeltaE1+DeltaE2)/3.0
axRelax.plot(sizes,-ave, 'ko', zorder=10)
axRelax.plot(sizes,DeltaE1-ave, 'ko', zorder=10)
axRelax.plot(sizes,DeltaE2-ave, 'ko', zorder=10)
axRelax.set_ylim(-4,4)


# ax.set_ylabel("$\Delta_{2}$, $\Delta_{3}$(meV)")
axRelax.set_xlabel("Size (nm)")
# ax.set_ylim(ymax = 20, ymin = 0)
plt.xlim(xmax = 10, xmin = 1)
axOrtho.set_title("Bright-Bright Splitting")
axOrtho.set_ylabel("Orthorhombic \n$\delta$ (meV)")
axRelax.set_ylabel("Relaxed \n$\delta$ (meV)")

fig.savefig("BrightBright.png", dpi=300, bbox_inches='tight')