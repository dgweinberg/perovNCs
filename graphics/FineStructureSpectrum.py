import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import sys


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


title = "Spectra"
print(title)

fig, (axCubic, axOrtho, axRelax) = plt.subplots(3,1,figsize = (5,7),sharex=True)


def plot_spectra(filename, ax):
    OS = np.loadtxt(filename)
    Ediff=OS[:,2]-OS[0,2]
    nstates = np.shape(OS)[0]
    xmin = -1
    xmax = 26
    xvals = np.linspace(xmin,xmax,1500)
    spectrum = np.zeros_like(xvals);

    spectrumx = np.zeros_like(xvals);
    spectrumy = np.zeros_like(xvals);
    spectrumz = np.zeros_like(xvals);
    sigma = 0.2
    nspect = np.shape(spectrum)[0]


    binnedOS = []
    binnedEne = []
    nDistStates = 0
    prevEne = -10

    for i in range(nstates):
        spectrum+=OS[i,3] *np.exp(-0.5*(xvals-Ediff[i]*27211.396641308)**2/(sigma**2))
        spectrumx+= 2.0/3.0 * OS[i,2]* (OS[i,4]**2+OS[i,5]**2) *np.exp(-0.5*(xvals-Ediff[i]*27211.396641308)**2/(sigma**2))
        spectrumy+= 2.0/3.0 * OS[i,2]* (OS[i,6]**2+OS[i,7]**2) *np.exp(-0.5*(xvals-Ediff[i]*27211.396641308)**2/(sigma**2))
        spectrumz+= 2.0/3.0 * OS[i,2]* (OS[i,8]**2+OS[i,9]**2) *np.exp(-0.5*(xvals-Ediff[i]*27211.396641308)**2/(sigma**2))
        if(Ediff[i]>(xmax+20)/27211.396641308):
            continue;
        if (Ediff[i]-prevEne > 0.1/27211.396641308):
                binnedOS.append(OS[i,3])
                binnedEne.append(Ediff[i]*27211.396641308)
                nDistStates+=1
                prevEne = Ediff[i]
        else:
                binnedOS[nDistStates-1]+=OS[i,3]
        
    binnedEne.append(0.00)
    binnedOS.append(0.00)
    binnedEne.append(1000.00)
    binnedOS.append(0.00)

    ax.set_yticks([])
    ax.set_xlim(xmin, xmax)
    ax.stem(binnedEne,np.array(binnedOS)+0.3*max(binnedOS), linefmt ='k--', markerfmt ='None', basefmt = 'k-', bottom=0.3*max(binnedOS), use_line_collection=True)
    ax.scatter(binnedEne,np.array(binnedOS)+0.3*max(binnedOS), c='Black', marker ='x', s=40)
    ax.plot(xvals,spectrum+0.3*max(binnedOS),'k-')
    ax.plot(xvals,spectrumx+0.2*max(binnedOS),'r-', label= 'x')
    ax.plot(xvals,spectrumy+0.1*max(binnedOS),'b-', label= 'y')
    ax.plot(xvals,spectrumz+0.0*max(binnedOS),'g-', label= 'z')
    


plot_spectra("../electronic/ortho/6x6x6/bse/OS.dat", axOrtho)
plot_spectra("../electronic/relax/6x6x6/bse/OS.dat", axRelax)
plot_spectra("../electronic/cubic/6x6x6/bse/OS.dat", axCubic)

# axOrtho.title(title)

#plt.legend()
fig.supxlabel("Energy (meV)", y=0.03)
fig.supylabel("Intensity (arb. units)",x=0.0)
fig.suptitle("Simulated Spectra", y=0.94)

axCubic.set_ylabel("Cubic")
axOrtho.set_ylabel("Orthorhombic")
axRelax.set_ylabel("Relaxed")
fig.savefig('FineStructure.png', dpi=300, bbox_inches='tight')


#plt.show()



