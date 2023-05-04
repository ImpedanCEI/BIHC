'''
This example shows how to compute the dissipated power
for an LHC measured spectrum retrieved from Timber

The impedance curve used is from a generic
Pillbox cavity simulated with CST studio
exported in .txt format

@date: Created on 03/05/2023
@author: Elena de la Fuente, Leonardo Sito
'''
import sys
sys.path.append('../../')

import bihc
import matplotlib.pyplot as plt
import numpy as np

try:
    import pytimber
except:
    print('This example uses pytimber. Please follow the installation guide to set it in your python environment')


def fillingSchemeLHC():
    '''Example filling scheme for LHC
    '''

    ninj = 10 # Defining number of injections
    nslots = 3564 # Defining total number of slots for LHC
    ntrain = 4 # Defining the number of trains
    nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 

    batchS = 7 # Batch spacing in 25 ns slots
    injspacing = 37 # Injection spacing in 25 ns slots

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchS
    stt = [False]*injspacing
    sc = [False]*(nslots-(ntrain*nbunches*ninj+((ntrain-1)*(batchS)*ninj)+((1)*injspacing*(ninj))))
    an1 = bt+ st +bt+ st+ bt+ st+ bt+ stt
    an = an1 * ninj + sc # This is the final true false sequence that is the beam distribution

    return an

# Set beam object with spectrum = 'user'

beam = bihc.Beam(Np=2.3e11, bunchLength=1.2e-9, fillingScheme=fillingSchemeLHC(), machine='LHC', spectrum='user') 

# Importing spectrum from timber
print('Accessing timber...')
db=pytimber.LoggingDB() 
time = pytimber.parsedate('2015-08-07 19:55:00')
data = db.get('ALB.SR4.B1:SPEC_POW', time)['ALB.SR4.B1:SPEC_POW'][1][0] #in dB
spectrum = 10.0**(data/20)
f = np.arange(0,len(spectrum))*beam.frev

#beam.setSpectrum([f, spectrum]) #???


plt.plot(spectrum)
plt.show()

'''
# Importing an impedance curve
impedance_file = 'PillboxImpedance.txt'

Z = bihc.Impedance()
Z.getImpedanceFromCST(impedance_file)

# built-in plot spectrum and normalized impedance
beam.plotSpectrumAndImpedance(Z)

# Computing the dissipated power value
print(f'Beam 6675 power loss: {beam.getPloss(Z)[0]} W')
'''