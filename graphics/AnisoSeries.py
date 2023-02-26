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

def plot_TT_Split(ax, all_sizes, fileDir,fileNames,fileExtension, fmtStr='k+', markerSize=1):
    sizes = []
    DeltaE1 = []
    DeltaE2 = []

    for n,file in enumerate(fileNames):
        print(f"Loading file: {fileDir+file+fileExtension}")
        try:
            OS = np.loadtxt(fileDir+file+fileExtension)
        except:
            print(f"Can't find: {fileDir+file+fileExtension}")
        else:
            if OS.size==0: continue;    
            OS[:,2]*=27211.396641308
            sizes.append(all_sizes[n])
            DeltaE1.append(OS[2,2]-OS[1,2])
            DeltaE2.append(OS[3,2]-OS[1,2])
    DeltaE1 = np.array(DeltaE1)
    DeltaE2 = np.array(DeltaE2)
    ave = (DeltaE1+DeltaE2)/3.0
    ax.plot(sizes,(DeltaE1+DeltaE2-ave)/3, fmtStr,alpha=1, markersize=markerSize)

    # ax.plot(sizes,-ave, fmtStr,alpha=1)
    # ax.plot(sizes,(DeltaE1-ave), fmtStr,alpha=1)
    # ax.plot(sizes,(DeltaE2-ave), fmtStr,alpha=1)
    return;

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

# my_4x4xN_sizes = Pb_to_nm(np.array([4*4*2, 
#                                     4*4*3,
#                                     4*4*4,
#                                     4*4*5,
#                                     4*4*6,
#                                     4*4*7,
#                                     4*4*8])**(1/3))

# my_5x5xN_sizes = Pb_to_nm(np.array([5*5*2, 
#                                     5*5*3,
#                                     5*5*4,
#                                     5*5*5,
#                                     5*5*6,
#                                     5*5*7,
#                                     5*5*8])**(1/3))


my_4x4xN_sizes = Pb_to_nm(np.array([4*4*3,
                                    4*4*4,
                                    4*4*5,
                                    4*4*6,
                                    4*4*7,
                                    4*4*8])**(1/3))

my_5x5xN_sizes = Pb_to_nm(np.array([5*5*3,
                                    5*5*4,
                                    5*5*5,
                                    5*5*6,
                                    5*5*7,
                                    5*5*8])**(1/3))

my_5x5xN_sizes = Pb_to_nm(np.array([6*6*3,
                                    6*6*4,
                                    6*6*5,
                                    6*6*6,
                                    6*6*7,
                                    6*6*8])**(1/3))

my_5x5xN_deg_of_aniso = np.array([3.,4.,5.,6.,7.,8.])/5.
my_4x4xN_deg_of_aniso = np.array([3.,4.,5.,6.,7.,8.])/4.
my_6x6xN_deg_of_aniso = np.array([3.,4.,5.,6.,7.,8.])/6.

orthoFileDir = "../electronic/ortho/"
relxFileDir = "../electronic//relax/"
cubeFileDir = "../electronic//cubic/"


fileExtension = "/bse/OS.dat"

my_5x5xN_names = ["5x5x3","4x5x5","5x5x5","5x5x6","5x5x7","5x5x8"]
my_4x4xN_names = ["3x4x4","4x4x4","4x4x5","4x4x6","4x4x7","4x4x8"]
my_6x6xN_names = ["6x6x3","6x6x4","5x6x6","6x6x6","6x6x7","6x6x8"]


fig, (ax1, ax2,ax3) = plt.subplots(3,1,figsize = (4,6),sharex=True)


plot_TT_Split(ax1,my_4x4xN_deg_of_aniso, cubeFileDir, my_4x4xN_names, fileExtension, 'bs',3)
plot_TT_Split(ax1,my_5x5xN_deg_of_aniso, cubeFileDir, my_5x5xN_names, fileExtension, 'bs',6)
plot_TT_Split(ax1,my_6x6xN_deg_of_aniso, cubeFileDir, my_6x6xN_names, fileExtension, 'bs',9)



plot_TT_Split(ax2,my_4x4xN_deg_of_aniso, orthoFileDir, my_4x4xN_names, fileExtension, 'gd',3)
plot_TT_Split(ax2,my_5x5xN_deg_of_aniso, orthoFileDir, my_5x5xN_names, fileExtension, 'gd',6)
plot_TT_Split(ax2,my_6x6xN_deg_of_aniso, orthoFileDir, my_6x6xN_names, fileExtension, 'gd',9)


plot_TT_Split(ax3,my_4x4xN_deg_of_aniso, relxFileDir, my_4x4xN_names, fileExtension, 'ko',3)
plot_TT_Split(ax3,my_5x5xN_deg_of_aniso, relxFileDir, my_5x5xN_names, fileExtension, 'ko',6)
plot_TT_Split(ax3,my_6x6xN_deg_of_aniso, relxFileDir, my_6x6xN_names, fileExtension, 'ko',9)





ax1.set_ylim([0,6])
ax2.set_ylim([0,6])
ax3.set_ylim([0,6])

ax2.set_ylabel("$\\bar{\Delta}$ (meV)")
ax3.set_xlabel("Aspect Ratio")
# ax.set_ylim(ymax = 20, ymin = 0)
# plt.xlim(xmax = 6, xmin = 1)
ax1.set_title("Bright-Bright Splitting")

fig.savefig("AnisoBrightBright.png", dpi=300, bbox_inches='tight')