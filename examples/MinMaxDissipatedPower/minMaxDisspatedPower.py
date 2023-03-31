'''
This example shows how to get the minimum and
maximum value of dissipated power when considering
a frequency uncertainty of +-20 Mhz

* date: 27/03/2023
* author: F. Giordano, E. de la Fuente, L. Sito
'''
import sys
sys.path.append('../../')

import bihc
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

#LPC csv file name
file = '25ns_2760b_2748_2494_2572_288bpi_13inj.csv'

# Data retrival from timber, with different bunch profile shapes
beam = bihc.Beam(LPCfile=file, bunchShape='GAUSSIAN', machine='LHC', Np=2.3e11, 
                 bunchLength=1.2e-9, fmax=2e9,spectrum='numeric', ppbk=250)
[f,S] = beam.spectrum

# Create Impedance object
impedance_file = 'PillboxImpedance.txt'
Z = bihc.Impedance()
Z.getImpedanceFromCST(impedance_file)

# Create 4 beam objects for each injection
Np = 1.2e11   # Number of protons per bunch
t0 = 25e-9    # Slot space [s]
bl = 1.65e-9  # total bunch length flat top

beam = bihc.Beam(bunchLength=bl, bunchShape='q-GAUSSIAN', qvalue=1.25, machine='SPS', fillMode='FLATTOP', fillingScheme=fillingSchemeSPS(4), Np=Np, d=t0) #4th injection flat top
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

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(beam.powerSpectrum[0], beam.powerSpectrum[1], 'b', label='Power spectrum')
ax.plot(fi, Zplot/np.max(Zplot), 'r', label='Shifted Impedance')
ax.plot(fi, Zi/np.max(Zi), 'k', label='Impedance')
ax.text(0.1, 0.85, f'Unshifted Power loss: {round(ploss0,2)} W', transform=ax.transAxes, color='k', weight='bold')
ax.text(0.1, 0.8, f'Max Shifted Power loss: {round(np.max(power),2)} W', transform=ax.transAxes, color='r', weight='bold')

ax.set_xlim(0, fmax)
ax.set_ylim(0, 1.0)
ax.set_xlabel('frequency [Hz]')
ax.set_ylabel('Amplitude [a.u.]')
ax.set_title('Maximum power loss with $f_{shift}$ of '+ f'{round(noise[np.argmax(power)]*deltaF/1e6,2)} MHz')
ax.legend()

fig.suptitle(f'SPS max power loss for Z: {impedance_file}')
plt.show()

