'''
Script to calculate BGI beam 
induced heating in the SPS

date: 30/03/23
author: edelafue, lsito, 
'''
import sys
sys.path.append('../../')

import bihc
import numpy as np
import matplotlib.pyplot as plt


def fillingSchemeSPS_standard(ntrains):
    '''
    Returns the filling scheme for the SPS

    Parameters
    ----------
    ntrains: number of trains (batches)
    nbunches: default, 48. number of bunches per train
    '''
    # Define filling scheme: parameters
    nslots = 920 # Defining total number of slots for SPS
    nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchspacing = 9 # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ntrains)
    an = (bt + st)*ntrains + sc

    return an


#select filling scheme
fillingSchemeSPS=fillingSchemeSPS_standard

# Reading Impedance file
Z = bihc.Impedance()
impedance_file = 'BGIfer1_wf_comp.txt'
#impedance_file = 'BGIv11_wf_re_im.txt'
Z.getImpedanceFromCST(impedance_file)
fmax = np.max(Z.f)

# Create 4 beam objects for each injection
Np = 1.2e11   # Number of protons per bunch
bl = 1.65e-9  # total bunch length flat top

b_bin = bihc.Beam(bunchLength=bl, bunchShape='BINOMIAL', exp=1.5, machine='SPS', fillMode='FLATTOP', fillingScheme=fillingSchemeSPS(4), Np=Np, fmax=fmax) #4th injection flat top
b_gauss = bihc.Beam(bunchLength=bl, bunchShape='GAUSSIAN', machine='SPS', fillMode='FLATTOP', fillingScheme=fillingSchemeSPS(4), Np=Np, fmax=fmax) #4th injection flat top
b_qgauss = bihc.Beam(bunchLength=bl, bunchShape='q-GAUSSIAN', qvalue=1.25, machine='SPS', fillMode='FLATTOP', fillingScheme=fillingSchemeSPS(4), Np=Np, fmax=fmax) #4th injection flat top

# get power
pwr_bin = b_bin.getPloss(Z)[0]
pwr_gauss = b_gauss.getPloss(Z)[0]
pwr_qgauss = b_qgauss.getPloss(Z)[0]

# Computing the dissipated power value
print(f'Computed power loss with binomial n=1.5: {pwr_bin} W')
print(f'Computed power loss with gaussian: {pwr_gauss} W')
print(f'Computed power loss with q-gaussian q=1.25: {pwr_qgauss} W')

fig, ax2 = plt.subplots(1,1, figsize=(8,6))

ax2.plot(b_gauss.powerSpectrum[0]/1e9, b_gauss.spectrum[1], c='b', label=f'Gaussian', alpha=0.7)
ax2.plot(b_qgauss.powerSpectrum[0]/1e9, b_qgauss.spectrum[1], label=f'q-Gaussiam q=1.5', c='g', alpha=0.7)
ax2.plot(b_bin.powerSpectrum[0]/1e9, b_bin.spectrum[1], label=f'Binomial n=1.25', c='r', alpha=0.7)

ax2.set_ylim(ymin=0)

axx2 = ax2.twinx()
axx2.plot(Z.f/1e9, Z.Zr, c='k', label='Impedance')
axx2.text(0.1, 0.9, f'Power loss Binomial: {round(pwr_bin,2)} W', transform=axx2.transAxes, color='k', weight='bold')
axx2.text(0.1, 0.85, f'Power loss Gaussian: {round(pwr_gauss,2)} W', transform=axx2.transAxes, color='k', weight='bold')
axx2.text(0.1, 0.8, f'Power loss q-Gaussian: {round(pwr_qgauss,2)} W', transform=axx2.transAxes, color='k', weight='bold')
axx2.set_ylabel('Impedance abs(Z) [$\Omega$]')
axx2.set_ylim(ymin=0)

ax2.set_xlim(0, fmax/1e9)
ax2.set_ylabel('Spectrum [a.u.]')
ax2.set_xlabel('frequency [GHz]')
ax2.legend()

fig.suptitle(f'SPS BGI power loss for Z: {impedance_file}')
fig.tight_layout()
plt.show()
