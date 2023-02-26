import matplotlib.pyplot as plt 
import numpy as np

#read expected band structure
exp = np.loadtxt('expBandStruct_0.dat');

#read best band structure
best = np.loadtxt('bandStruct_0.dat')


fig, ax = plt.subplots();

nbands=np.shape(exp)[1]
for i in range(0,nbands-1):
	if(i==1):
		ax.plot(exp[:,0], exp[:,1+i], 'r.', label='Literature');
		ax.plot(best[:,0], best[:,1+i], 'k-', label='Best Fit');
	else:
		ax.plot(exp[:,0], exp[:,1+i], 'r.', label='');
		ax.plot(best[:,0], best[:,1+i], 'k-',  label='');

ax.plot(best[:11,0],-5.92638 - 27.2114/(2.0*0.24)*best[:11,0]**2, '-m',  label='Effective Mass')
ax.plot(best[:11,0],-4.15521 + 27.2114/(2.0*0.24)*best[:11,0]**2, '-m',  label='')

ax.plot(best[30:,0],-5.92638 - 27.2114/(2.0*0.24)*(best[30:,0]-best[40,0])**2, '-m',  label='')
ax.plot(best[30:,0],-4.15521 + 27.2114/(2.0*0.24)*(best[30:,0]-best[40,0])**2, '-m',  label='')

ax.set_title('Orthorhombic Band Structure')
ax.set_xlabel('K-Point')
ax.set_ylabel('Energy (eV)')
##uncomment this to set up labels on the special symmetry points
##set up the labels at the special points
#need to know where the high symmetry points are in the k-list
gammaIndex=0
yIndex=10
tIndex=20
zIndex=30
gamma2Index=40
xIndex=50
#re-write this with the points you go through in the k-path
plt.xticks(ticks=[exp[gammaIndex,0], exp[yIndex,0], exp[tIndex,0], exp[zIndex,0], exp[gamma2Index,0], exp[xIndex,0]], labels=[r'$\Gamma$', 'Y','T', 'Z',r'$\Gamma$', 'X']);

plt.xlim(exp[0,0],exp[-1,0])
plt.ylim(-8,-2)
plt.legend()
plt.savefig('Plot.png',format='png');

