'''
This example compares different bunch profile shapes 
using a same filling scheme given by an LPC csv file and 
the same impedance curve.

It plots the impact of the different bunch shapes in the 
beam spectrum and computes the difference in power loss

* date: 12/12/2022
* author: Francesco Giordano, Elena de la Fuente, Leonardo Sito
'''

import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
import bihc

#LPC beam filling scheme 
file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'

# Data retrival from timber, with different bunch profile shapes
b_gauss= bihc.Beam(LPCfile=file, bunchShape='GAUSSIAN', ppbk=1000, verbose=False)
b_qgauss= bihc.Beam(LPCfile=file, bunchShape='q-GAUSSIAN', ppbk=1000, verbose=False)
b_bin = bihc.Beam(LPCfile=file, bunchShape='BINOMIAL', ppbk=1000,  verbose=False)
b_cos = bihc.Beam(LPCfile=file, bunchShape='COS2', ppbk=1000, verbose=False)
b_par = bihc.Beam(LPCfile=file, bunchShape='PARABOLIC', ppbk=1000, verbose=False)

# Storing spectra 
[f_gauss, S_gauss] = b_gauss.spectrum
[f_qgauss, S_qgauss] = b_qgauss.spectrum
[f_bin, S_bin] = b_bin.spectrum
[f_cos, S_cos] = b_cos.spectrum
[f_par, S_par] = b_par.spectrum

# Storing profile
[t_gauss, s_gauss] = b_gauss.profile_1_bunch
[t_qgauss, s_qgauss] = b_qgauss.profile_1_bunch
[t_bin, s_bin] = b_bin.profile_1_bunch
[t_cos, s_cos] = b_cos.profile_1_bunch
[t_par, s_par] = b_par.profile_1_bunch

# Plotting 
fig, ax = plt.subplots(1, figsize=(12,7))

ax.plot(t_gauss, s_gauss, 'b', lw=3, label='gaussian')
ax.plot(t_qgauss, s_qgauss, 'g', lw=3, label='q-gaussian')
ax.plot(t_cos, s_cos, c='m', lw=3, label='cosine')
ax.plot(t_bin, s_bin,c='orange', lw=3, label='binomial')
ax.plot(t_par, s_par, c='red', lw=3, label='parabolic')
ax.legend()

fig, axs = plt.subplots(5,1, figsize=(8,12), dpi=200)

axs[0].plot(f_gauss, S_gauss,'b')
axs[1].plot(f_qgauss, S_qgauss,'g')
axs[2].plot(f_bin, S_bin, c='orange' )
axs[3].plot(f_cos, S_cos, c='m')
axs[4].plot(f_par, S_par, c='red')

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
power_qgauss = b_qgauss.getPloss(Z)[0]
power_bin = b_bin.getPloss(Z)[0]
power_cos = b_cos.getPloss(Z)[0]
power_par = b_par.getPloss(Z)[0]

axs[0].text(0.3, 0.8, f'Gaussian Power loss: {round(power_gauss,2)} W', transform=axs[0].transAxes, color='k', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
axs[1].text(0.3, 0.8, f'q-Gaussian Power loss: {round(power_qgauss,2)} W', transform=axs[1].transAxes, color='k', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
axs[2].text(0.3, 0.8, f'Binomial Power loss: {round(power_bin,2)} W', transform=axs[2].transAxes, color='k', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
axs[3].text(0.3, 0.8, f'Cosine Squared Power loss: {round(power_cos,2)} W', transform=axs[3].transAxes,color='k', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
axs[4].text(0.3, 0.8, f'Parabolic Power loss: {round(power_par,2)} W', transform=axs[4].transAxes,color='k', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))

print(f'Gaussian power loss: {power_gauss} W')
print(f'Gaussian power loss: {power_qgauss} W')
print(f'Binomial power loss: {power_bin} W')
print(f'Cosine Squared power loss: {power_cos} W')
print(f'Parabolic power loss: {power_par} W')

plt.tight_layout()
plt.show()
