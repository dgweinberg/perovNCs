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
                                    9*9*9, 10*10*10, 11*11*11])**(1/3))

my_cubic_size = Pb_to_nm(np.array([2*2*2, 
                                    3*3*3,
                                    4*4*4,
                                    5*5*5,
                                    6*6*6,
                                    7*7*7,
                                    8*8*8,
                                    9*9*9, 10*10*10, 11*11*11])**(1/3))


my_rect_size = Pb_to_nm(np.array([2*2*3, 2*3*3, 
                                     3*3*4, 3*4*4,
                                     4*4*5, 4*5*5, 
                                     5*5*6, 5*6*6, 
                                     6*6*7, 6*7*7,
                                     7*7*8, 7*8*8])**(1/3))


cubeFileDir = "../electronic/relax/"
rectFileDir = "../electronic/production2/relax/"
spinfileExtension = "/bse/spins.dat"
excitonfileExtension = "/bse/exciton.dat"
OSfileExtension = "/bse/OS.dat"

cubeDirNames = ["2x2x2", "3x3x3", "4x4x4", "5x5x5", "6x6x6", "7x7x7", "8x8x8", "9x9x9"]#, "10x10x10", "11x11x11"]
rectDirNames = ["2x2x3", "2x3x3", "3x3x4",
                "3x4x4", "4x4x5", "4x5x5", "5x5x6", "5x6x6",
                "6x6x7", "6x7x7", "7x7x8", "7x8x8"]

fig, (ax,ax2) = plt.subplots(1,2,figsize = (9,4))
fig.tight_layout(pad=2.0)

for n,file in enumerate(cubeDirNames):
    print(f"Loading files from: {cubeFileDir+file}")
    try:
        spins = np.loadtxt(cubeFileDir+file+spinfileExtension, usecols=range(0,5))
        exciton = np.loadtxt(cubeFileDir+file+excitonfileExtension, skiprows=1)
        OS = np.loadtxt(cubeFileDir+file+OSfileExtension)
    except:
        print(f"Can't find all files in: {cubeFileDir+file}")
    else: 

    	if(n==0): 
    		ax.scatter(my_cubic_size[n]*np.ones_like(spins[:4,4]),spins[:4,4],c = "Black")
    		ax.scatter(my_cubic_size[n]*np.ones_like(spins[4:7,4]),spins[4:7,4],c = "Red")
    		ax2.scatter(my_cubic_size[n]*np.ones_like(spins[:4,4]),np.divide(27.2114*exciton[:4,4],exciton[:4,6]),c ="Black")
    		ax2.scatter(my_cubic_size[n]*np.ones_like(spins[4:7,4]),np.divide(27.2114*exciton[4:7,4],exciton[4:7,6]),c ="Red")
    		#ax2.scatter(my_cubic_size[n]*np.ones_like(spins[:7,4]),np.divide(-27.2114*exciton[:7,3],exciton[:7,6]))
    	else: 
    		ax.scatter(my_cubic_size[n]*np.ones_like(spins[0,4]),spins[0,4],c ="Black")
    		ax.scatter(my_cubic_size[n]*np.ones_like(spins[1:4,4]),spins[1:4,4],c ="Red")
    		ax2.scatter(my_cubic_size[n]*np.ones_like(spins[0,4]),np.divide(27.2114*exciton[0,4],exciton[0,6]),c ="Black")
    		ax2.scatter(my_cubic_size[n]*np.ones_like(spins[1:4,4]),np.divide(27.2114*exciton[1:4,4],exciton[1:4,6]),c ="Red")

ax.set_title("Exciton Spin")
ax.set_ylabel(r"$\langle \hat{S}_{\textrm{tot}}^2 \rangle_n$")
ax.set_xlabel("Size (nm)")
ax.set_ylim(1.25,2)
ax.set_xlim(1,8)

ax2.set_title("Exchange Interaction")
ax2.set_ylabel(r"$\frac{\langle \hat{K}_{x} \rangle_n}{E_B^n}$")
ax2.set_xlabel("Size (nm)")
ax2.set_xlim(1,8)
# plt.show()


#plt.show()
fig.savefig("ExcitonStats.png", dpi=300, bbox_inches='tight')
