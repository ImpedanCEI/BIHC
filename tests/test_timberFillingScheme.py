'''
This example shows how to compute the dissipated power
for an LHC fillNumber 

The filling scheme and the bunch parameters are 
obtained from a fill number through Timber.
The impedance curve used is from a generic
Pillbox cavity simulated with CST studio
exported in .txt format

@date: Created on 08/12/2022
@author: Elena de la Fuente, Leonardo Sito
'''
import sys
sys.path.append('../../')

import bihc

# Data retrival from timber, with different bunch profile shapes
b_6675 = bihc.Beam(fillNumber=6675, bunchShape='GAUSSIAN')
# Exporting frequency array, needed for the impedance curve
[f, S] = b_6675.spectrum

# Importing an impedance curve
impedance_file = 'PillboxImpedance.txt'

Z = bihc.Impedance()
Z.getImpedanceFromCST(impedance_file)

# built-in plot spectrum and normalized impedance
b_6675.plotSpectrumAndImpedance(Z)

# Computing the dissipated power value
print(f'Beam 6675 power loss: {b_6675.getPloss(Z)[0]} W')

