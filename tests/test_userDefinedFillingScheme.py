'''
This example shows the definition of a filling 
scheme as a list of True or Falses values.
It is then used along some parameters of the bunches 
(Number of protons and slot space) to define a Beam object.
Using some built in methods, the time distribution 
and beam spectrum are plotted. The manual plotting 
is also included below.
Finally, the power loss is calculated for a simple Pillbox 
impedance map

* date: 12/12/2022
* author: Francesco Giordano, Elena de la Fuente, Leonardo Sito
'''
import sys
sys.path.append('../../')

import bihc
import matplotlib.pyplot as plt

# Define filling scheme: parameters
ninj = 10 # Defining number of injections
nslots = 3564 # Defining total number of slots for LHC
ntrain = 4 # Defining the number of trains
nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 

batchS = 7 # Batch spacing in 25 ns slots
injspacing = 37 # Injection spacing in 25 ns slots
BS = 200 # Batch spacing in ns

Np = 1.2e11 # Number of protons per bunch
t0 = 25e-9 # Slot space [s]

# Defining the trains as lists of True/Falses
bt = [True]*nbunches
st = [False]*batchS
stt = [False]*injspacing
sc = [False]*(nslots-(ntrain*nbunches*ninj+((ntrain-1)*(batchS)*ninj)+((1)*injspacing*(ninj))))
an1 = bt+ st +bt+ st+ bt+ st+ bt+ stt
an = an1 * ninj + sc # This is the final true false sequence that is the beam distribution

# Data retrival from timber
custom_beam = bihc.Beam(bunchShape='GAUSSIAN', beamNumber=1, fillingScheme=an, Np=Np, d=t0, verbose=False)

# built-in plotting
custom_beam.plotLongitudinalProfile()
custom_beam.plotPowerSpectrum()

# Retrieve attributes
[t,profile] = custom_beam.longitudinalProfile
[f,spectrum] = custom_beam.spectrum
[f,pspectrum] = custom_beam.powerSpectrum

# Manual plotting
fig, (ax1,ax2) = plt.subplots(2,1)

ax1.plot(t*1e9, profile, c='b', label='profile')
ax1.set_ylabel('Intensity [p/b]')
ax1.set_xlabel('time [ns]')
ax1.legend()

ax2.plot(f/1e9, spectrum, c='b', label='spectrum')
ax2.plot(f/1e9, pspectrum, c='r', label='power spectrum')
ax2.set_xlim(0, 2)
ax2.set_ylabel('normalized amplitude')
ax2.set_xlabel('frequency [GHz]')
ax2.legend()

fig.suptitle('User defined filling scheme')
fig.set_size_inches(12,6)
fig.tight_layout()
plt.show()

# Adding power loss
Z = bihc.Impedance(f)
Z.getImpedanceFromCST('PillboxImpedance.txt')

# Computing the dissipated power value
print(f'Custom beam power loss: {custom_beam.getPloss(Z)[0]} W')