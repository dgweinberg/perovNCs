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



par = np.array([0.0, 0.25, 0.50, 0.75, 1.0])


nPPVBM = np.array([-6.21681, -6.12131, -6.07166, -6.06753, -6.10502])
nPPCBM = np.array([-4.28850, -4.39066, -4.42299, -4.39210, -4.31613])

expectedVBM = np.array([-6.25, np.nan, np.nan, np.nan, -6.13])
expectedCBM = np.array([-4.47, np.nan, np.nan, np.nan, -4.41])

DFTVBM = np.array([3.3598, 3.4787, 3.5283, 3.5192, 3.4676])
DFTCBM = np.array([4.6932, 4.7975, 4.9050, 5.0113, 5.1098])

VBMass = np.array([-0.374785, -0.360157, -0.347122, -0.335590, -0.325492])
CBMass = np.array([0.424322,  0.413935,  0.404260,  0.395348,  0.387253])


matplotlib.rcParams.update({'figure.figsize' : (6,10)})


fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
plt.setp((ax1, ax2, ax3), xticks = [0,0.25,0.5,0.75,1], xticklabels=['cubic', '', '', '', 'orthorhombic'])



ax1.set_title("DFT Band Edges")
ax1.plot(par,DFTCBM, 'bo-', label ="CBM")
ax1.plot(par,DFTVBM, 'bo-', label ="VBM")
ax1.set_ylabel("Energy (ev)")
# ax1.legend()

ax2.set_title("Pseudopotential Band Edges")

ax2.plot(par,nPPCBM, 'ro-', label ="CBM")
ax2.plot(par,nPPVBM, 'ro-', label ="VBM")
ax2.plot(par, expectedCBM, 'ko', label = "Expected")
ax2.plot(par, expectedVBM, 'ko')
ax2.set_ylabel("Energy (ev)")


ax3.set_title("Band Gaps")
ax3.plot(par, DFTCBM-DFTVBM, 'bo-', label = "DFT")
ax3.plot(par, nPPCBM-nPPVBM, 'ro-', label = "Pseudopotential")
ax3.plot(par,expectedCBM-expectedVBM, 'ko', label= "Expected")
ax3.set_ylabel("Energy (ev)")





plt.savefig("mixBandsCompDFT.png",dpi=300, bbox_inches='tight')
