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

xsize = 5;
ysize = xsize*3/4

orthoI_s = np.loadtxt("../electronic/pots/ortho/perovBandsOrtho-35245/final_run/pot_I.dat")
orthoI_j32 = np.loadtxt("../electronic/pots/ortho/perovBandsOrtho-35245/final_run/pot_I_j3_2.dat")
orthoI_j12 = np.loadtxt("../electronic/pots/ortho/perovBandsOrtho-35245/final_run/pot_I_j1_2.dat")

cubicI_s = np.loadtxt("../electronic/pots/cubic/perovBandsCubic-34960-09-22/finalRun/pot_I.dat")
cubicI_j32 = np.loadtxt("../electronic/pots/cubic/perovBandsCubic-34960-09-22/finalRun/pot_I_j3_2.dat")
cubicI_j12 = np.loadtxt("../electronic/pots/cubic/perovBandsCubic-34960-09-22/finalRun/pot_I_j1_2.dat")

fig, ax = plt.subplots(figsize = (xsize,ysize))

ax.plot(orthoI_s[:,0], orthoI_s[:,1], 'k--')
ax.plot(orthoI_j12[:,0], orthoI_j12[:,1], 'g--')
ax.plot(orthoI_j32[:,0], orthoI_j32[:,1], 'r--')



ax.plot(cubicI_s[:,0], cubicI_s[:,1], 'k:')
ax.plot(cubicI_j32[:,0], cubicI_j32[:,1], 'r:')
ax.plot(cubicI_j12[:,0], cubicI_j12[:,1], 'g:')


ax.plot([0],[0], 'k-', label="$s$")
ax.plot([0],[0], 'g-', label="$p_{\\frac{1}{2}}$")
ax.plot([0],[0], 'r-', label="$p_{\\frac{3}{2}}$")
plt.legend()

plt.ylabel("Energy (Hartree)")
plt.xlabel("Radius (Bohr)")
# plt.ylim(ymax = 3.5, ymin = 1.5)
plt.xlim(xmax = 5, xmin = 0)
plt.title("Iodine Pseudopotentials")

#plt.show()
plt.savefig("IPots.png", dpi=300, bbox_inches='tight')

######################################################################################################

orthoI_s = np.loadtxt("../electronic/pots/ortho/perovBandsOrtho-35245/final_run/pot_Pb.dat")
orthoI_j32 = np.loadtxt("../electronic/pots/ortho/perovBandsOrtho-35245/final_run/pot_Pb_j3_2.dat")
orthoI_j12 = np.loadtxt("../electronic/pots/ortho/perovBandsOrtho-35245/final_run/pot_Pb_j1_2.dat")

cubicI_s = np.loadtxt("../electronic/pots/cubic/perovBandsCubic-34960-09-22/finalRun/pot_Pb.dat")
cubicI_j32 = np.loadtxt("../electronic/pots/cubic/perovBandsCubic-34960-09-22/finalRun/pot_Pb_j3_2.dat")
cubicI_j12 = np.loadtxt("../electronic/pots/cubic/perovBandsCubic-34960-09-22/finalRun/pot_Pb_j1_2.dat")

fig, ax = plt.subplots(figsize = (xsize,ysize))



ax.plot(orthoI_s[:,0], orthoI_s[:,1], 'k--')
ax.plot(orthoI_j32[:,0], orthoI_j32[:,1], 'r--')
ax.plot(orthoI_j12[:,0], orthoI_j12[:,1], 'g--')


ax.plot(cubicI_s[:,0], cubicI_s[:,1], 'k:')
ax.plot(cubicI_j32[:,0], cubicI_j32[:,1], 'r:')
ax.plot(cubicI_j12[:,0], cubicI_j12[:,1], 'g:')

ax.plot([0],[0], 'k-', label="$s$")
ax.plot([0],[0], 'g-', label="$p_{\\frac{1}{2}}$")
ax.plot([0],[0], 'r-', label="$p_{\\frac{3}{2}}$")
plt.legend()

plt.ylabel("Energy (Hartree)")
plt.xlabel("Radius (Bohr)")
# plt.ylim(ymax = 3.5, ymin = 1.5)
plt.xlim(xmax = 5, xmin = 0)
plt.title("Lead Pseudopotentials")
#plt.show()
plt.savefig("PbPots.png", dpi=300, bbox_inches='tight')