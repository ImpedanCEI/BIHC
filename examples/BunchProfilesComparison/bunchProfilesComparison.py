'''
This example compares different bunch profile shapes 
using a same filling scheme given by an LPC csv file and 
the same impedance curve.

It plots the impact of the different bunch shapes in the 
beam spectrum and power loss result

* date: 12/12/2022. Revised 31/01/2024
* author: Elena de la Fuente, Leonardo Sito
'''

import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
import bihc

#LPC beam filling scheme 
file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'

# Create beam object
bl = 1.2e-9                 # bunch length [s]
Np = 2.3e11                 # bunch intensity [protons/bunch]
fillMode = 'FLATTOP'        # Energy
fmax = 2e9                  # Maximum frequency of the beam spectrum [Hz]
ppbk = 250 					# number of samples per slot
verbose = False 				# Enable terminal verbosy output 

b_gauss = bihc.Beam(Np=Np, bunchLength=bl, LPCfile=file, bunchShape='GAUSSIAN', ppbk=ppbk, verbose=verbose)
b_qgauss = bihc.Beam(Np=Np, bunchLength=bl, LPCfile=file, bunchShape='q-GAUSSIAN', qvalue=1.2, ppbk=ppbk, verbose=verbose)
b_bin = bihc.Beam(Np=Np, bunchLength=bl, LPCfile=file, bunchShape='BINOMIAL', exp=2.5, ppbk=ppbk, verbose=verbose)
b_cos = bihc.Beam(Np=Np, bunchLength=bl, LPCfile=file, bunchShape='COS2', ppbk=ppbk, verbose=verbose)
b_par = bihc.Beam(Np=Np, bunchLength=bl, LPCfile=file, bunchShape='PARABOLIC', ppbk=ppbk, verbose=verbose)

#  ------- Plotting in time ------------

# Store 1 bunch profile (t)
[t_gauss, s_gauss] = b_gauss.profile_1_bunch
[t_qgauss, s_qgauss] = b_qgauss.profile_1_bunch
[t_bin, s_bin] = b_bin.profile_1_bunch
[t_cos, s_cos] = b_cos.profile_1_bunch
[t_par, s_par] = b_par.profile_1_bunch

fig, ax = plt.subplots(1, figsize=(8,6), dpi=150)
ax.plot(t_gauss, s_gauss, 'b', lw=3, label='gaussian')
ax.plot(t_qgauss, s_qgauss, 'g', lw=3, label='q-gaussian')
ax.plot(t_cos, s_cos, c='m', lw=3, label='cosine')
ax.plot(t_bin, s_bin,c='orange', lw=3, label='binomial')
ax.plot(t_par, s_par, c='red', lw=3, label='parabolic')
ax.legend()
fig.suptitle('Bunch shape comparison in time for 1 bunch slot')
fig.tight_layout()
plt.show()

# Integral check: Area under profile == 1
dt = t_gauss[2]-t_gauss[1]
print('\nNormalization check: time integral == 1.0')
print(f'Gaussian integral: {np.sum(s_gauss)*dt} ')
print(f'q-Gaussian integral: {np.sum(s_qgauss)*dt} ')
print(f'Binomial integral: {np.sum(s_bin)*dt} ')
print(f'Cosine Squared integral: {np.sum(s_cos)*dt} ')
print(f'Parabolic integral: {np.sum(s_par)*dt} ')

#   Importing an impedance curve (Resonator)
# ----------------------------------------------
Z = bihc.Impedance()
Z.getResonatorImpedance(Rs=7e3, Qr=1e2, fr=1.75e9, f=b_gauss.spectrum[0])

# ---------  Plotting in frequency -------------
# Storing spectra 
[f_gauss, S_gauss] = b_gauss.spectrum
[f_qgauss, S_qgauss] = b_qgauss.spectrum
[f_bin, S_bin] = b_bin.spectrum
[f_cos, S_cos] = b_cos.spectrum
[f_par, S_par] = b_par.spectrum

fig, axs = plt.subplots(5,1, figsize=(8,10), dpi=100)
for ax in axs:
	ax.plot(Z.f, Z.Zr/np.max(Z.Zr), 'k')

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

print('\n- Bunch shape impact in specturm and power for a narrowband impedance' )
print(f'Gaussian power loss for resonator impedance: {power_gauss} W')
print(f'q-Gaussian power loss for resonator impedance: {power_qgauss} W')
print(f'Binomial power loss for resonator impedance: {power_bin} W')
print(f'Cosine Squared power loss for resonator impedance: {power_cos} W')
print(f'Parabolic power loss for resonator impedance: {power_par} W')

fig.suptitle('Bunch shape impact in specturm and power for a narrowband impedance' , fontweight='bold')
fig.tight_layout()
plt.show()

# Importing an impedance curve (Resistive wall)
# ----------------------------------------------
Z = bihc.Impedance()
Z.getRWImpedance(L=1.0, b=15e-3, sigma=5.7e7, f=f_gauss) 

# ---------  Plotting in frequency -------------
fig, axs = plt.subplots(5,1, figsize=(8,10), dpi=100)

for ax in axs:
	ax.plot(Z.f, Z.Zr/np.max(Z.Zr), 'k')

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

print('\n- Bunch shape impact in specturm and power for a broadband impedance' )
print(f'Gaussian power loss for resistive wall impedance: {power_gauss} W')
print(f'q-Gaussian power loss for resistive wall impedance: {power_qgauss} W')
print(f'Binomial power loss for resistive wall impedance: {power_bin} W')
print(f'Cosine Squared power loss for resistive wall impedance: {power_cos} W')
print(f'Parabolic power loss for resistive wall impedance: {power_par} W')

fig.suptitle('Bunch shape impact in specturm and power for a broadband impedance', fontweight='bold')
fig.tight_layout()
plt.show()

