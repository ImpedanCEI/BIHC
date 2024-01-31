'''
This example shows how to get the minimum and
maximum value of dissipated power when considering
a frequency uncertainty of +-20 Mhz

* date: 27/03/2023. Revised: 31/01/2024
* author: E. de la Fuente, L. Sito
'''
import sys
sys.path.append('../../')

import bihc
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

# SPS user defined filling scheme
def fillingSchemeSPS_standard(ntrains):
    '''
    Returns the filling scheme for the SPS

    Parameters
    ----------
    ntrains: number of injections (batches)
    '''
    # Define filling scheme: parameters
    ntrain = 1 # SPS has 1 train per cycle
    nslots = 924 # Defining total number of slots for SPS
    nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchspacing = 9 # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ntrains)
    an = (bt + st)*ntrains + sc

    return an


#------- Ploss calculation ----------

# Create beam object
fillingScheme = fillingSchemeSPS_standard(ntrains=4)
bl = 1.2e-9                 # bunch length [s]
Np = 2.3e11                 # bunch intensity [protons/bunch]
bunchShape = 'q-GAUSSIAN'     # bunch profile shape in time 
qvalue = 1.2                # value of q parameter in the q-gaussian distribution
fillMode = 'FLATTOP'        # Energy
fmax = 2e9                  # Maximum frequency of the beam spectrum [Hz]

beam = bihc.Beam(Np=Np, bunchLength=bl, fillingScheme=fillingScheme,
                bunchShape=bunchShape, qvalue=qvalue, 
                machine='SPS', fillMode=fillMode, spectrum='numeric', fmax=fmax)
[f,S] = beam.spectrum

# Create Impedance object
impedance_file = 'PillboxImpedance.txt'
Z = bihc.Impedance(f)
Z.getImpedanceFromCST(impedance_file)

# Get unshifted ploss 
ploss, ploss_density = beam.getPloss(Z) 

#---------------- Rigid Shift power loss ------------------------------
shift = 40e6  # distance between shift steps [Hz]
shifts, power = beam.getShiftedPloss(Z, shift=shift)

print(f'Minimum dissipated power: P_min = {np.min(power)}, at step {shifts[np.argmin(power)]}')
print(f'Maximum dissipated power: P_max = {np.max(power)}, at step {shifts[np.argmax(power)]}')
print(f'Average dissipated power: P_mean = {np.mean(power)}')

#-----------   Plotting  ---------
Zmax = beam.Zmax                        # Shifted impedance object
f = beam.powerSpectrum[0]               # frequency array
pow_spectrum = beam.powerSpectrum[1]    # Power spectrum Λ^2

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(f, pow_spectrum, 'g', label='Power spectrum')

axx = ax.twinx()
axx.plot(Zmax.f, Zmax.Zr/np.max(Zmax.Zr), 'r', label='Shifted Impedance')
axx.plot(Z.f, Z.Zr/np.max(Z.Zr), 'k', label='Impedance')
axx.set_ylabel('Normalized impedance [a.u]', color='r')
axx.set_ylim(0,1.2)
axx.legend()

ax.text(0.05, 0.95, f'Unshifted Power loss: {round(ploss,2)} W', transform=ax.transAxes, color='k', weight='bold')
ax.text(0.05, 0.9, f'Max Shifted Power loss: {round(np.max(power),2)} W', transform=ax.transAxes, color='r', weight='bold')

ax.set_xlim(0, fmax)
ax.set_ylim(0, 1.2)
ax.set_xlabel('frequency [Hz]')
ax.set_ylabel('Power spectrum amplitude [a.u.]')
ax.set_title('Maximum power loss with $f_{shift}$ of '+ f'{round(shifts[np.argmax(power)]*shift/1e6,2)} MHz')

fig.suptitle(f'SPS power loss for Z="{impedance_file}"')
plt.show()

