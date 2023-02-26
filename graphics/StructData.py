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

my_sizes = Pb_to_nm(np.array([2*2*2, 
                                    3*3*3,
                                    4*4*4,
                                    5*5*5,
                                    6*6*6,
                                    7*7*7,
                                    8*8*8,
                                    9*9*9, 10*10*10, 11*11*11])**(1/3))

cubic_bond_angles = np.array([180, 179.995729, 180, 179.997865, 180, 
 179.998576, 179.99878, 180, 179.999051, 179.999146])

cubic_DLC = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

relax_bond_angles = np.array([168.921413, 171.87659, 167.880618, 166.967248, 165.446422,
 164.172722, 163.155145, 162.286901 ,161.554343, 160.199934])

relax_DLC = np.array([0.00444, 0.01955, 0.09757, 0.12845, 0.12414, 0.088933, 0.059757, 0.04555, 0.021156, 0.02298])

ortho_bond_angles = np.array([154.107144, 154.107222, 154.107238, 154.107305, 154.107077,
 154.107121, 154.107162, 154.107119, 154.10715, 154.107122])

ortho_DLC = np.array([0.02614, 0.02614, 0.02614, 0.02615, 0.02616, 0.026167, 0.026157, 0.02615, 0.026133, 0.02614])

DFT_bond_angles = np.array([165.054, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
DFT_DLC = np.array([0.000017, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

Bulk_bond_angles = np.array([156.101])
Bulk_DLC = np.array([0.0286])

fig, ax = plt.subplots(figsize = (4,3))

ax.scatter(my_sizes, cubic_bond_angles, c='Blue', marker='s', label='Cubic')
ax.scatter(my_sizes, DFT_bond_angles, c='Red', marker='x', label='DFT', s=80)
ax.scatter(my_sizes, ortho_bond_angles, c='Green', marker='d', label='Orthorhombic')
ax.scatter(7.5,Bulk_bond_angles,c='Black', marker='x', label='Bulk', s=80)
ax.scatter(my_sizes, relax_bond_angles, c='Black', marker='o', label='Relaxed')


ax.set_ylabel("Pb-I-Pb Bond Angle (Deg.)")
ax.set_xlabel("Size (nm)")
ax.set_ylim(ymax = 185, ymin = 143)
plt.xlim(xmax = 8, xmin = 1)
# ax.set_title("Relaxed Structures")
#plt.legend(fontsize='small', frameon=False, loc=(-0.05,0.0), ncols=3, columnspacing=0.0, handletextpad=0.0)
fig.savefig("BondAngles.png", dpi=300, bbox_inches='tight')

################################################################

fig, ax = plt.subplots(figsize = (4,3))

ax.scatter(my_sizes, cubic_DLC, c='Blue', marker='s', label='Cubic')

ax.scatter(my_sizes, ortho_DLC, c='Green', marker='d', label='Orthorhombic')

ax.scatter(my_sizes, relax_DLC, c='Black', marker='o', label='Relaxed')

ax.scatter(my_sizes, DFT_DLC, c='Red', marker='x', label='DFT',s=80)

ax.scatter(7.5,Bulk_DLC,c='Black', marker='x', label='Bulk', s=80)

ax.set_ylabel(r"Lattice Anisotropy (\AA)")
ax.set_xlabel("Size (nm)")
ax.set_ylim(ymin = -0.04)
plt.xlim(xmax = 8, xmin = 1)
plt.legend(fontsize='small', frameon=False, loc=(1,0.3), ncols=1, columnspacing=0.0, handletextpad=0.0)
fig.savefig("LatticeAnisotropy.png", dpi=300, bbox_inches='tight')



