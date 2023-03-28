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
Zreal=np.interp(f,Z.f,Z.Zr)
Z.f = f

#Initialice variables
deltaF = Z.f[1] - Z.f[0]
power = np.array([])
steps = np.array([], dtype=int)

# Create noise for +-20 MHz uncertainty
size = 2000
noise = 20e6 * np.random.uniform(low=-1,size=200)  #random uniform (Francesco & Niccolo)
noise = 20e6 * np.linspace(-1,1,size)              #controlled step

# Start scan
print(f'Starting scan with {size} points...')
for i, n in tqdm(enumerate(noise), "Computing scan: ", total=size):
    step = int(n/deltaF)
    steps = np.append(steps, step)
    Z.Zr = np.roll(Zreal, step) 
    Z.Zr[0] = 0.0
    power = np.append(power, beam.getPloss(Z)[0])

print(f'Minimum dissipated power: P_min = {np.min(power)}, at step {np.argmin(power)}')
print(f'Maximum dissipated power: P_max = {np.max(power)}, at step {np.argmax(power)}')
print(f'Average dissipated power: P_mean = {np.mean(power)}')

# plotting
step = steps[np.argmax(power)]
Z.Zr = np.roll(Zreal, step) 
Z.Zr[:step] = Z.Zr[:step]*0.0

fig, ax = plt.subplots()
ax.plot(beam.powerSpectrum[0], beam.powerSpectrum[1], 'b', label='Power spectrum')
ax.plot(Z.f, Z.Zr/np.max(Z.Zr), 'r', label='shifted Impedance')
ax.plot(f, Zreal/np.max(Zreal), 'k', label='Impedance')

ax.set_xlim(0, 2e9)
ax.set_ylim(0,1.0)
ax.set_xlabel('frequency [Hz]')
ax.set_ylabel('Amplitude [a.u.]')
ax.set_title('Maximum power loss with $f_{shift}$ of '+ f'{round(noise[np.argmax(power)]/1e6,2)} MHz')
ax.legend()
plt.show()

