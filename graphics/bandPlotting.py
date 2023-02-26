import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
fig, ax = plt.subplots(figsize = (xsize,ysize))



#read expected band structure
exp = np.loadtxt('../electronic/pots/ortho/perovBandsOrtho-35245/final_run/expBandStruct_0.dat');

#read best band structure
best = np.loadtxt('../electronic/pots/ortho/perovBandsOrtho-35245/final_run/bandStruct_0.dat')


nbands=np.shape(exp)[1]
for i in range(0,nbands-1):
	if(i==1):
		ax.plot(exp[:,0], exp[:,1+i], 'r.', label='Literature');
		ax.plot(best[:,0], best[:,1+i], 'k-', label='Best Fit');
	else:
		ax.plot(exp[:,0], exp[:,1+i], 'r.', label='');
		ax.plot(best[:,0], best[:,1+i], 'k-',  label='');

# ax.plot(best[:11,0],-6.25 - 27.2114/(2.0*0.24)*best[:11,0]**2, '-m',  label='Effective Mass')
# ax.plot(best[:11,0],-4.32 + 27.2114/(2.0*0.24)*best[:11,0]**2, '-m',  label='')


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
plt.ylim(-7,-3)
# plt.legend()
plt.savefig('OrthoBands.png', dpi=300, bbox_inches='tight');

##############################################################################################################################################

fig, ax = plt.subplots(figsize = (xsize,ysize))



#read expected band structure
exp = np.loadtxt('../electronic/pots/cubic/perovBandsCubic-34960-09-22/finalRun/expBandStruct_0.dat');

#read best band structure
best = np.loadtxt('../electronic/pots/cubic/perovBandsCubic-34960-09-22/finalRun/bandStruct_0.dat')


nbands=np.shape(exp)[1]
for i in range(0,nbands-1):
	if(i==1):
		ax.plot(exp[:,0], exp[:,1+i], 'r.', label='Literature');
		ax.plot(best[:,0], best[:,1+i], 'k-', label='Best Fit');
	else:
		ax.plot(exp[:,0], exp[:,1+i], 'r.', label='');
		ax.plot(best[:,0], best[:,1+i], 'k-',  label='');

# ax.plot(best[:11,0],-6.25 - 27.2114/(2.0*0.24)*best[:11,0]**2, '-m',  label='Effective Mass')
# ax.plot(best[:11,0],-4.32 + 27.2114/(2.0*0.24)*best[:11,0]**2, '-m',  label='')


ax.set_title('Cubic Band Structure')
ax.set_xlabel('K-Point')
ax.set_ylabel('Energy (eV)')
#uncomment this to set up labels on the special symmetry points
#set up the labels at the special points
#need to know where the high symmetry points are in the k-list
rIndex=0
gammaIndex=11
xIndex=21
mIndex=31
#re-write this with the points you go through in the k-path
plt.xticks(ticks=[exp[rIndex,0], exp[gammaIndex,0], exp[xIndex,0], exp[mIndex,0]], labels=['R', r'$\Gamma$','X', 'M']);

plt.xlim(exp[0,0],exp[-1,0])
plt.ylim(-8,-1)
# plt.legend()
plt.savefig('CubicBands.png', dpi=300, bbox_inches='tight');