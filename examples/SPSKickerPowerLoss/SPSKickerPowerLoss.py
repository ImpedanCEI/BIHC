'''
This example shows a power loss computation
for the MKPL kicker in the CERN SPS for a beam
with 4 trains and 48 bunches 25ns spaced and 
1.5e11 intensity p/b. 

Based on the studies by C. Zannini.

* date: 01/03/2023
* author: F. Giordano, E. de la Fuente, L. Sito
'''
import sys
sys.path.append('../../')

import bihc
import numpy as np
import matplotlib.pyplot as plt

def fillingSchemeSPS(ninj, nubuches=48):
    '''
    Returns the filling scheme for the SPS

    Parameters
    ----------
    ninj: number of injections (batches)
    nbunches: default, 48. number of bunches per train
    '''
    # Define filling scheme: parameters
    ntrain = 1 # SPS has 1 train per cycle
    nslots = 920 # Defining total number of slots for SPS
    nbunches = 48 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchspacing = 8 # Batch spacing in 25 ns slots (200 ns)

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ninj)
    an = (bt + st)*ninj + sc

    return an

# Reading Impedance file
Z = bihc.Impedance()
Z.getImpedanceFromCST('ZserHVreference.txt')
fmax = np.max(Z.f)

# Create 4 beam objects for each injection
Np = 1.5e11   # Number of protons per bunch
t0 = 25e-9    # Slot space [s]
bl = 7.505192141958421e-10 * 4 #0.225 m / clight *4 (to fit Francesco's definition: -2sigma, +2sigma)

b1 = bihc.Beam(bunchLength=bl, machine='SPS', fillMode='FB', fillingScheme=fillingSchemeSPS(1), Np=Np, d=t0, fmax=fmax) #first injection
b2 = bihc.Beam(bunchLength=bl, machine='SPS', fillMode='FB', fillingScheme=fillingSchemeSPS(2), Np=Np, d=t0, fmax=fmax) #second injection
b3 = bihc.Beam(bunchLength=bl, machine='SPS', fillMode='FB',  fillingScheme=fillingSchemeSPS(3), Np=Np, d=t0, fmax=fmax) #third injection
b4 = bihc.Beam(bunchLength=bl, machine='SPS', fillMode='FB', fillingScheme=fillingSchemeSPS(4), Np=Np, d=t0, fmax=fmax) #4th injection flat bottom
b4ft = bihc.Beam(bunchLength=bl, machine='SPS', fillMode='FLATTOP', fillingScheme=fillingSchemeSPS(4), Np=Np, d=t0, fmax=fmax) #4th injection flat top

# Computing the dissipated power value
avg_powerloss = (b1.getPloss(Z)[0] + b2.getPloss(Z)[0] + b3.getPloss(Z)[0] + b4.getPloss(Z)[0] + b4ft.getPloss(Z)[0])/5
print(f'Computed power loss: {avg_powerloss} W')

# Plotting
fig, (axs) = plt.subplots(4,1, figsize=(4,6))

axs[3].plot(b4.longitudinalProfile[0]*1e9, b4.longitudinalProfile[1], label='4 trains')
axs[2].plot(b3.longitudinalProfile[0]*1e9, b3.longitudinalProfile[1], label='3 trains')
axs[1].plot(b2.longitudinalProfile[0]*1e9, b2.longitudinalProfile[1], label='2 trains')
axs[0].plot(b1.longitudinalProfile[0]*1e9, b1.longitudinalProfile[1], label='1 train')

for ax in axs:
    ax.set_xlim(0, t0*920*1e9)
    ax.set_ylabel('Intensity [p/b]')
    ax.set_xlabel('time [ns]')
    ax.legend()

fig.suptitle('SPS bunches considered')
fig.tight_layout()
plt.show()

fig, ax2 = plt.subplots(1,1, figsize=(6,4))

ax2.plot(b4.powerSpectrum[0]/1e9, b4.spectrum[1], label='4 trains')
ax2.plot(b3.powerSpectrum[0]/1e9, b3.spectrum[1], label='3 trains')
ax2.plot(b2.powerSpectrum[0]/1e9, b2.spectrum[1], label='2 trains')
ax2.plot(b1.powerSpectrum[0]/1e9, b1.spectrum[1], label='1 train')
ax2.set_ylim(ymin=0)

axx2 = ax2.twinx()
axx2.plot(Z.f/1e9, abs(Z.Z), c='k', label='Impedance')
#axx2.fill_between(Z.f/1e9, abs(Z.Z), color='k', alpha=0.2)
axx2.text(0.5, 0.2, f'Max Power loss: {round(b4ft.getPloss(Z)[0],2)} W', transform=axx2.transAxes, color='k', weight='bold')
axx2.text(0.5, 0.1, f'Avg Power loss: {round(avg_powerloss,2)} W', transform=axx2.transAxes, color='k', weight='bold')
axx2.set_ylabel('Impedance abs(Z) [$\Omega$]')
axx2.set_ylim(0, 3000)

ax2.set_xlim(0, fmax/1e9)
ax2.set_ylabel('Spectrum [a.u.]')
ax2.set_xlabel('frequency [GHz]')
ax2.legend()

fig.suptitle('SPS Kicker power loss')
fig.tight_layout()
plt.show()
