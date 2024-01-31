'''
Example to read filling schemes generated
by the LPC tool, in .csv format

It plots the comparison of three different
filling shcemes and computes the power loss

@date: Created on 08/12/2022. Revised 31/01/24.
@author: Leonardo Sito, Elena de la Fuente
'''

import sys
sys.path.append('../../') #if bihc is not pip installed

import matplotlib.pyplot as plt
import bihc

# LPC csv file names
# downloaded from LPC: https://lpc.web.cern.ch/cgi-bin/fillingSchemeTab.py
file1 = '25ns_2748b_2736_2258_2374_288bpi_12inj.csv'
file2 = '25ns_2374b_2361_1730_1773_236bpi_13inj_hybrid_2INDIV.csv' 
file3 = '8b4e_1972b_1967_1178_1886_224bpi_12inj.csv'

#------- Ploss calculation ----------

# Create beam object
bl = 1.2e-9  		# bunch length [s]
Np = 2.3e11			# bunch intensity [protons/bunch]
bunchShape = 'GAUSSIAN' 	# bunch profile shape in time 
fillMode = 'FLATTOP' 		# Energy
fmax = 2e-9					# Maximum frequency of the beam spectrum [Hz]

beamBcms1 = bihc.Beam(LPCfile=file1, Np=Np, bunchLength=bl, bunchShape=bunchShape, fillMode=fillMode, fmax=fmax)
beamBcms2 = bihc.Beam(LPCfile=file2,  Np=Np, bunchLength=bl, bunchShape=bunchShape, fillMode=fillMode,fmax=fmax)
beam8b4e = bihc.Beam(LPCfile=file3,  Np=Np, bunchLength=bl,  bunchShape=bunchShape, fillMode=fillMode, fmax=fmax)

# Importing an impedance curve
impedance_file = 'PillboxImpedance.txt'
Z = bihc.Impedance()
Z.getImpedanceFromCST(impedance_file)

# Computing the dissipated power value for the different filling schemes
print(f'25ns_2748b_2736_2258_2374_288bpi_12inj power loss: {beamBcms1.getPloss(Z)[0]} W')
print(f'25ns_2760b_2748_2494_2572_288bpi_13inj power loss: {beamBcms2.getPloss(Z)[0]} W')
print(f'8b4e_1972b_1967_1178_1886_224bpi_12inj power loss: {beam8b4e.getPloss(Z)[0]} W')

#------- Plotting ----------

# Storing spectra 
[f_bcms1, S_bcms1] = beamBcms1.spectrum
[f_bcms2, S_bcms2] = beamBcms2.spectrum
[f_8b4e, S_8b4e] = beam8b4e.spectrum

fig, axs = plt.subplots(3,1, figsize=[8,10])

# Plotting beam spectrum
axs[0].plot(f_bcms1, S_bcms1, label=beamBcms1.LPCfile, c='b')
axs[1].plot(f_bcms2, S_bcms2, label=beamBcms2.LPCfile, c='g')
axs[2].plot(f_8b4e, S_8b4e,label=beam8b4e.LPCfile, c='darkorange')

axs[0].set_title(label=beamBcms1.LPCfile, color='b')
axs[1].set_title(label=beamBcms2.LPCfile, color='g')
axs[2].set_title(label=beam8b4e.LPCfile, color='darkorange')

for ax in axs:
	# Plotting impedance
	axx = ax.twinx()
	axx.plot(Z.f, Z.Zr, c='r')

	axx.set_ylabel('Impedance Re(Z) [$\Omega$]', color='r')
	ax.set_xlim(0, 2e9)
	ax.set_xlabel('frequency [Hz]')
	ax.set_ylabel('Spectrum [a.u.]')

plt.tight_layout()
plt.show()