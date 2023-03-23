'''
Example to read filling schemes generated
by the LPC tool, in .csv format

It plots the comparison of three different
filling shcemes and computes the power loss

@date: Created on 08/12/2022
@author: Leonardo Sito, Elena de la Fuente
'''

import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import bihc

#LPC csv file names
file1='25ns_2748b_2736_2258_2374_288bpi_12inj.csv'
file2='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
file3='8b4e_1972b_1967_1178_1886_224bpi_12inj.csv'

# Data retrival from timber, with different bunch profile shapes
beamBcms1 = bihc.Beam(LPCfile=file1, bunchShape='GAUSSIAN')
beamBcms2 = bihc.Beam(LPCfile=file2, bunchShape='GAUSSIAN')
beam8b4e = bihc.Beam(LPCfile=file3, bunchShape='GAUSSIAN')

# Storing spectra 
[f_bcms1, S_bcms1] = beamBcms1.spectrum
[f_bcms2, S_bcms2] = beamBcms2.spectrum
[f_8b4e, S_8b4e] = beam8b4e.spectrum

# Plotting 
fig, axs = plt.subplots(3,1)

axs[0].plot(f_bcms1, S_bcms1, label=beamBcms1.LPCfile)
axs[1].plot(f_bcms2, S_bcms2, label=beamBcms2.LPCfile)
axs[2].plot(f_8b4e, S_8b4e,label=beam8b4e.LPCfile)

for ax in axs:
	ax.set_xlim(0, 2e9)
	ax.set_xlabel('frequency [Hz]')
	ax.set_ylabel('Spectrum [a.u.]')
	ax.legend()

plt.tight_layout()
plt.show()

# Importing an impedance curve
impedance_file = 'PillboxImpedance.txt'
Z = bihc.Impedance()
Z.getImpedanceFromCST(impedance_file)

# Computing the dissipated power value for the different Bunch Profiles
print(f'25ns_2748b_2736_2258_2374_288bpi_12inj power loss: {beamBcms1.getPloss(Z)[0]} W')
print(f'25ns_2760b_2748_2494_2572_288bpi_13inj power loss: {beamBcms2.getPloss(Z)[0]} W')
print(f'8b4e_1972b_1967_1178_1886_224bpi_12inj power loss: {beam8b4e.getPloss(Z)[0]} W')

