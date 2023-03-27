'''
This example shows how to get the minimum and
maximum value of dissipated power when considering
a frequency uncertainty of +-20 Mhz

* date: 27/03/2023
* author: F. Giordano, E. de la Fuente, L. Sito
'''
import sys
sys.path.append('../../')

import bihc
import numpy as np
import matplotlib.pyplot as plt


#LPC csv file names
file = '25ns_2748b_2736_2258_2374_288bpi_12inj'

# Data retrival from timber, with different bunch profile shapes
beam = bihc.Beam(LPCfile=file1, bunchShape='GAUSSIAN') #TODO add rest of important parameters

# Create Impedance object
impedance_file = 'PillboxImpedance.txt'
Z = Impedance(f)
Z.getImpedanceFromCST(impedance_file)

#Initialice variables
deltaF = Z.f[1] - Z.f[0]
Power = []
Power_density = []

# Create noise for +-20 MHz uncertainty
size = 1000
noise = 20e6 * np.random.uniform(low=-1,size=size) #random uniform
noise = 40e6 * np.linspace(-1,1,size)             #controlled step

# Start scan
Print(f'Starting scan with {size} points...')
for i, n in enumerate(noise):

    step = int(n/deltaF)

    Z.Zr = np.roll(Z.real, step) #TODO is not roll

    Power.append(beam.getPloss(Z))
    Power_density.append(beam.getPloss(Z))

Print(f'Maximum dissipated power: P_max = {np.max(Power)}')
Print(f'Average dissipated power: P_mean = {np.mean(Power)}')