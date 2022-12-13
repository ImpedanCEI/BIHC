'''
In this example there is the definition of a filling scheme from a fill number.
The filling scheme and the parameters of the bunch are retrived from Timber.
The only parameter that still needs to be set by the user is the bunch profile
shape. The main differences among different bunch shapes is shown through a plot.

@date: Created on 08/12/2022
@author: lsito
'''
# Example: using an LHC fill number to obtain beam parameters
import sys
sys.path.append('../../../')

import bihc

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Data retrival from timber, with different bunch profile shapes
b_6675_gauss = bihc.Beam(fillNumber=6675, bunchShape='GAUSSIAN')
b_6675_bin = bihc.Beam(fillNumber=6675, bunchShape='BINOMIAL')
b_6675_cos = bihc.Beam(fillNumber=6675, bunchShape='COS2')

# Storing spectra 
[f_gauss, S_gauss] = b_6675_gauss.spectrum
[f_bin, S_bin] = b_6675_bin.spectrum
[f_cos, S_cos] = b_6675_cos.spectrum

# Plotting 
fig, axs = plt.subplots(3,1)

axs[0].plot(f_gauss, S_gauss)
axs[1].plot(f_bin, S_bin)
axs[2].plot(f_cos, S_cos)

plt.tight_layout()
plt.show()

# Importing an impedance curve
impedance_file = 'Impedance_file.txt'

Z = bihc.Impedance(f_gauss)
Z.getImpedanceFromCST(impedance_file)

# Computing the dissipated power value for the different Bunch Profiles
print(f'Gaussian power loss: {b_6675_gauss.getPloss(Z)[0]} W')
print(f'Binomial power loss: {b_6675_bin.getPloss(Z)[0]} W')
print(f'Cosine Squared power loss: {b_6675_cos.getPloss(Z)[0]} W')

