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


#https://pubs.acs.org/doi/full/10.1021/acs.nanolett.1c03114
son_size = np.array([6.5])
son_gap = np.array([1.849])


#https://aip.scitation.org/doi/full/10.1063/1.5109894
lian_size = np.array([11.8,10.8,6.5])
lian_errorbars = np.array([1.6,1.5,1.1])
lian_gap = 1239.84193/ np.array([676.0,665.0,650.0]) 


#https://pubs.acs.org/doi/full/10.1021/jacs.8b11447
yao_size = np.array([14.06, 11.79, 10.1, 7.33, 5.38, 4.70])
yao_gap = np.array([1.75, 1.76, 1.77, 1.79, 1.80, 1.83])

#https://pubs.acs.org/doi/full/10.1021/acs.jpcc.9b05969
demir_size = np.array([15])
demir_errorbars = np.array([3])
demir_gap = 1239.84193/np.array([683])

#https://www.science.org/doi/10.1126/science.aag2700
luther_size = np.array([3.4, 4.5, 5.0, 6.8, 8.0, 9.0, 12.5])
luter_gap =  1239.84193/ np.array([598, 616, 636, 648, 661, 668, 681])

#https://www.sciencedirect.com/science/article/pii/S092534672100848X
wang_size = np.array([36])
wang_gap = np.array([1.79])


#https://pubs.acs.org/doi/full/10.1021/acs.jpclett.2c01463
samanta_size = np.array([9.0,12.5])
samanta_errorbars = np.array([1.0,2.0])
samanta_gap = 1239.84193/np.array([677,685])

#https://www.nature.com/articles/nature25147
efros_size = np.array([14,14,14,14])
efros_gap = np.array([1.827,1.832,1.778,1.785])




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


cubeFileDir = "/home/dweinberg/Projects/rashba/calcs/production2/cubic/"
rectFileDir = "/home/dweinberg/Projects/rashba/calcs/production2/cubic/"
fileExtension = "/bse/OS.dat"

cubeDirNames = ["2x2x2", "3x3x3", "4x4x4", "5x5x5", "6x6x6", "7x7x7", "8x8x8", "9x9x9", "10x10x10", "11x11x11"]
rectDirNames = ["2x2x3", "2x3x3", "3x3x4",
                "3x4x4", "4x4x5", "4x5x5", "5x5x6", "5x6x6",
                "6x6x7", "6x7x7", "7x7x8", "7x8x8"]

fig, ax = plt.subplots(figsize = (4,3))

for n,file in enumerate(cubeDirNames):
    print(f"Loading file: {cubeFileDir+file+fileExtension}")
    try:
        OS = np.loadtxt(cubeFileDir+file+fileExtension)
    except:
        print(f"Can't find: {cubeFileDir+file+fileExtension}")
    else:
        if (np.size(OS)==0): continue    
        OS[:,2]*=27.211396641308

        ax.scatter(my_cubic_size[n]*np.ones_like(OS[0,2]), OS[0,2], c = "Black", marker = 's', alpha=1)
        ax.scatter(my_cubic_size[n]*np.ones_like(OS[1:3,2]), OS[1:3,2], c ='Red', marker = 's', alpha=1)

# ###################################
ax.errorbar(son_size,son_gap,fmt ='.', color='grey',label="Son")
ax.errorbar(lian_size,lian_gap,xerr=lian_errorbars,fmt = '.',color='grey', label="Lian")
ax.errorbar(yao_size,yao_gap, fmt = '.', color='grey',label="Yao")
ax.errorbar(demir_size, demir_gap, xerr=demir_errorbars, fmt = '.',color='grey', label = "Demir")
ax.errorbar(luther_size,luter_gap, fmt = '.', color='grey',label = "Luther")
ax.errorbar(wang_size,wang_gap, fmt = '.', color='grey',label = "Wang")
ax.errorbar(samanta_size,samanta_gap, xerr=samanta_errorbars, fmt= '.',color='grey', label = "Samanta")
ax.errorbar(efros_size, efros_gap, fmt = '.',color='grey', label="Efros")
# ###################################






plt.ylabel("Excitation Energy (eV)")
plt.xlabel("Size (nm)")
plt.ylim(ymax = 3.5, ymin = 1.5)
plt.xlim(xmax = 10, xmin = 1)
plt.title("Cubic Structures")
#plt.show()
plt.savefig("PerfectOS.png", dpi=300, bbox_inches='tight')
