'''
Script to calculate BGI beam 
induced heating in the SPS.

This calculates the maximum and 
minimum power loss possible for 
a q-gaussian beam

date: 30/03/23
author: edelafue, lsito
'''
import sys
sys.path.append('../../')

import bihc
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm


def fillingSchemeSPS_standard(ntrains):
    '''
    Returns the filling scheme for the SPS

    Parameters
    ----------
    ntrains: number of injections (batches)
    '''
    # Define filling scheme: parameters
    ntrain = 1 # SPS has 1 train per cycle
    nslots = 924 # Defining total number of slots for SPS
    nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchspacing = 9 # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ntrains)
    an = (bt + st)*ntrains + sc

    return an

def fillingSchemeSPS_8b4e(ntrains):
    '''
    Returns the filling scheme for the SPS 
    using the 8b4e pattern

    Parameters
    ----------
    ntrains: number of trains (batches)
    '''
    # Define filling scheme: parameters
    nslots = 924 # Defining total number of slots for SPS
    nbunches = 8*7 # Defining number of bunches e.g. 18, 36, 72.. 
    nempty = 4*6  #Defining number of empty slots between bunches
    batchspacing = 8 # Batch spacing in 25 ns slots (200 ns)

    # Defining the trains as lists of True/Falses
    bt = ([True]*8+[False]*4)*6+[True]*8
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+nempty+batchspacing)*ntrains)
    an = (bt + st)*ntrains + sc

    return an

#select filling scheme
fillingSchemeSPS=fillingSchemeSPS_standard

# Reading Impedance file
Z = bihc.Impedance()
impedance_file = 'impedance_with_ferrite.txt'
impedance_file = 'SPSWiresWFerritesSide.txt'
Z.getImpedanceFromCST(impedance_file)

# Create 4 beam objects for each injection
Np = 1.8e11   # Number of protons per bunch
bl = 1.55e-9  # total bunch length flat top

beam18 = bihc.Beam(bunchLength=bl, Np=1.8e11, bunchShape='GAUSSIAN', machine='SPS', fillMode='FLATTOP', fillingScheme=fillingSchemeSPS(4)) #4th injection flat top
beam23 = bihc.Beam(bunchLength=bl, Np=2.3e11, bunchShape='GAUSSIAN', machine='SPS', fillMode='FLATTOP', fillingScheme=fillingSchemeSPS(4)) #4th injection flat top

beam = beam23
[f,S] = beam.spectrum
ploss0 = beam.getPloss(Z)[0]

#---------------- Max power loss ------------------------------

#Initialice variables
size = 1000  #number of frequency samples in a 20MHz range
fmax = np.minimum(Z.f.max(), f.max())
fi = np.linspace(0, fmax, int(size*fmax/20e6))
deltaF = fi[1] - fi[0]

#interpolate
Zi = np.interp(fi, Z.f, Z.Zr) #interpolate to match beam frequencies
Zi[Zi<0] = 0.0
Z.f, Z.Zr = fi, Zi

# Create noise for +-20 MHz uncertainty
noise = np.arange(-size, size, 1, dtype=int)  #controlled step
power = np.array([])

# Start scan
print(f'Starting scan with {2*size} points, {np.round(deltaF/1e3,2)} kHz freq. shift...')

for step in tqdm(noise, "Computing scan: ", total=2*size):
    Z.Zr = np.roll(Zi, step) 
    if step > 0: Z.Zr[:step] = 0.0
    if step < 0: Z.Zr[:2*size-step] = 0.0
    power = np.append(power, beam.getPloss(Z)[0])

print(f'Minimum dissipated power: P_min = {np.min(power)}, at step {noise[np.argmin(power)]}')
print(f'Maximum dissipated power: P_max = {np.max(power)}, at step {noise[np.argmax(power)]}')
print(f'Average dissipated power: P_mean = {np.mean(power)}')

# plotting
step = noise[np.argmax(power)]
Zplot= np.roll(Zi, step) 
if step > 0: Zplot[:step] = 0.0
if step < 0: Zplot[:2*size-step] = 0.0

fig, ax = plt.subplots(figsize=(8,6), dpi=150)
ax.plot(beam.powerSpectrum[0], beam.powerSpectrum[1], 'b', label='Power spectrum')
ax.text(0.2, 0.85, f'Unshifted Power loss: {round(ploss0,2)} W', transform=ax.transAxes, color='k', weight='bold')
ax.text(0.2, 0.9, f'Max Shifted Power loss: {round(np.max(power),2)} W', transform=ax.transAxes, color='r', weight='bold')

axx = ax.twinx()
axx.plot(fi, Zplot, 'r', label='Shifted Impedance')
axx.plot(fi, Zi, 'k', label='Impedance')
axx.set_ylim(ymin=0)
axx.set_ylabel('Impedance')

ax.set_xlim(0, fmax)
ax.set_ylim(0, 1.0)
ax.set_xlabel('frequency [Hz]')
ax.set_ylabel('Spectrum mplitude [a.u.]')
ax.set_title('Maximum power loss with $f_{shift}$ of '+ f'{round(noise[np.argmax(power)]*deltaF/1e6,2)} MHz', fontsize=10)
#ax.legend()

fig.suptitle(f'SPS WS power loss with ferrite \n Intensity = 2.3e11 p/b, Bunch length = 1.55 ns \n')
plt.show()
