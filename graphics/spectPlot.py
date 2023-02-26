import numpy as np
import math
import matplotlib.pyplot as plt
import sys

title = sys.argv[1]
print(title)
OS = np.loadtxt('OS.dat')

nstates = np.shape(OS)[0]
xmin = np.around(OS[0,2]*27211.396641308, -2) - 100
xmax = 30000
xvals = np.linspace(xmin,xmax,15000)
spectrum = np.zeros_like(xvals);

spectrumx = np.zeros_like(xvals);
spectrumy = np.zeros_like(xvals);
spectrumz = np.zeros_like(xvals);
sigma = 20.0
nspect = np.shape(spectrum)[0]


binnedOS = []
binnedEne = []
nDistStates = 0
prevEne = -10

for i in range(nstates):
    spectrum+=OS[i,3] *np.exp(-0.5*(xvals-OS[i,2]*27211.396641308)**2/(sigma**2))
    spectrumx+= 2.0/3.0 * OS[i,2]* (OS[i,4]**2+OS[i,5]**2) *np.exp(-0.5*(xvals-OS[i,2]*27211.396641308)**2/(sigma**2))
    spectrumy+= 2.0/3.0 * OS[i,2]* (OS[i,6]**2+OS[i,7]**2) *np.exp(-0.5*(xvals-OS[i,2]*27211.396641308)**2/(sigma**2))
    spectrumz+= 2.0/3.0 * OS[i,2]* (OS[i,8]**2+OS[i,9]**2) *np.exp(-0.5*(xvals-OS[i,2]*27211.396641308)**2/(sigma**2))
    if(OS[i,2]>3.0/27.211396641308):
        continue;
    if (OS[i,2]-prevEne > 1.0/27211.396641308):
            binnedOS.append(OS[i,3])
            binnedEne.append(OS[i,2]*27211.396641308)
            nDistStates+=1
            prevEne = OS[i,2]
    else:
            binnedOS[nDistStates-1]+=OS[i,3]
    
binnedEne.append(0.00)
binnedOS.append(0.00)
fig = plt.gcf()
fig, ax = plt.subplots(figsize=(6,3))
ax.set_xlim(xmin, 500+xmin)
ax.stem(binnedEne,binnedOS, linefmt ='k--', markerfmt ='kx', basefmt = 'k-', use_line_collection=True)
ax.plot(xvals,spectrum,'k-')
ax.plot(xvals,spectrumx,'r--', label= 'x')
ax.plot(xvals,spectrumy,'b-.', label= 'y')
ax.plot(xvals,spectrumz,'g:', label= 'z')
plt.title(title)

plt.legend()


ax.set_xlabel("Energy (meV)")
ax.set_ylabel("Intensity (arb. units)")
fig.savefig('spectrum.png', dpi=300, bbox_inches='tight')


#plt.show()



