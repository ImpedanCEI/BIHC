'''
This example compares different bunch profile shapes 
using a same filling scheme given by an LPC csv file and 
the same impedance curve.

It plots the impact of the different bunch shapes in the 
beam spectrum and computes the difference in power loss

@date: Created on 08/12/2022
@author: Elena de la Fuente, Leonardo Sito
'''

import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import bihc

#LPC beam filling scheme 
file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'

# Data retrival from timber, with different bunch profile shapes
b_gauss= bihc.Beam(LPCfile=file, bunchShape='GAUSSIAN')
b_bin = bihc.Beam(LPCfile=file, bunchShape='BINOMIAL')
b_cos = bihc.Beam(LPCfile=file, bunchShape='COS2')

# Storing spectra 
[f_gauss, S_gauss] = b_gauss.spectrum
[f_bin, S_bin] = b_bin.spectrum
[f_cos, S_cos] = b_cos.spectrum

# Plotting 
fig, axs = plt.subplots(3,1)

axs[0].plot(f_gauss, S_gauss)
axs[1].plot(f_bin, S_bin)
axs[2].plot(f_cos, S_cos)

for ax in axs:
	ax.set_xlim(0, 2e9)
	ax.set_ylim(0, 1)
	ax.set_xlabel('frequency [Hz]')
	ax.set_ylabel('Spectrum [a.u.]')

# Importing an impedance curve
impedance_file = 'PillboxImpedance.txt'
Z = bihc.Impedance(f_gauss)
Z.getImpedanceFromCST(impedance_file)

# Computing the dissipated power value for the different Bunch Profiles
power_gauss = b_gauss.getPloss(Z)[0]
power_bin = b_bin.getPloss(Z)[0]
power_cos = b_cos.getPloss(Z)[0]

axs[0].text(0.3, 0.8, f'Gaussian Power loss: {round(power_gauss,2)} W', transform=axs[0].transAxes, color='tab:blue', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
axs[1].text(0.3, 0.8, f'Binomial Power loss: {round(power_bin,2)} W', transform=axs[1].transAxes, color='tab:blue', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
axs[2].text(0.3, 0.8, f'Cosine Squared Power loss: {round(power_cos,2)} W', transform=axs[2].transAxes,color='tab:blue', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))

print(f'Gaussian power loss: {power_gauss} W')
print(f'Binomial power loss: {power_bin} W')
print(f'Cosine Squared power loss: {power_cos} W')

plt.tight_layout()
plt.show()
