'''
This example shows how to compute the dissipated power. The filling scheme and
the bunch parameters are obtained from a fill number through Timber.
The impedance curve of the accelerator component is imported from a txt file
exported from CST.

@date: Created on 08/12/2022
@author: lsito
'''
# Example: using an LHC fill number to obtain beam parameters
import bihc

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Data retrival from timber, with different bunch profile shapes
b_6675 = bihc.Beam(fillNumber=6675, bunchShape='GAUSSIAN')
# Exporting frequency array, needed for the impedance curve
[f, S] = b_6675.spectrum

# Importing an impedance curve
impedance_file = 'Impedance_file.txt'

Z = bihc.Impedance(f)
Z.getImpedanceFromCST(impedance_file)

# Computing the dissipated power value
print(f'Beam 6675 power loss: {b_6675.getPloss(Z)[0]} W')

