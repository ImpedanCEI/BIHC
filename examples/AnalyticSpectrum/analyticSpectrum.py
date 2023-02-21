'''
This example compares the dissipated power difference obtained between 
using a numeric spactrum calculation or the analytic formula by C.Zannini
using a same filling scheme given by an `LPC Tool` generated `.csv` file and 
the same impedance curve from a simple Pillbox resonator.

It plots the impact of the spectrum type chosen in the 
beam spectrum and computes the difference in power loss

* date: 21/02/2022
* author: Francesco Giordano, Elena de la Fuente, Leonardo Sito
'''

import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
import bihc


# Beam data with different bunch profile shapes from LPC beam filling scheme 
file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
profiles = ['GAUSSIAN'] #, 'PARABOLIC', 'COS2'] #only GAUSSIAN works so far
power={}

# Plotting 
fig, axs = plt.subplots(len(profiles),1)

for i, prof in enumerate(profiles):

    beam_numeric = bihc.Beam(LPCfile=file, bunchShape=prof, verbose=False, spectrum='numeric')
    beam_analytic  = bihc.Beam(LPCfile=file, bunchShape=prof, verbose=False, spectrum='analytic')

    # Storing spectra 
    [fn, Sn] = beam_numeric.spectrum
    [fa, Sa] = beam_analytic.spectrum

    # Storing profile
    [tn, sn] = beam_numeric.profile_1_bunch
    [ta, sa] = beam_analytic.profile_1_bunch

    # plot spectrum
    if len(profiles) > 1: ax = axs[i]
    else: ax = axs

    ax.plot(fa, Sa, 'r+-', label='analytic')
    ax.plot(fn, Sn, 'bo-', label='numeric')

    # plot spectrum envelope
    sa_i = beam_analytic.lambdas[1]
    ax.plot(fa, sa_i, 'r', alpha=0.6)

    # Compute power loss

    # Importing an impedance curve
    impedance_file = 'PillboxImpedance.txt'
    Z = bihc.Impedance(fn)
    Z.getImpedanceFromCST(impedance_file)

    powern = beam_numeric.getPloss(Z)[0]
    powera = beam_analytic.getPloss(Z)[0]

    ax.text(0.3, 0.8, f'Numeric Power loss: {round(powern,2)} W', transform=ax.transAxes, color='tab:blue', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
    ax.text(0.3, 0.7, f'Analytic Power loss: {round(powera,2)} W', transform=ax.transAxes, color='tab:red', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))

    power[prof] = {'numeric': powern, 'analytic': powera}

    ax.legend()
    ax.set_xlim(0, 2e9)
    ax.set_ylim(0, 1)
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('Spectrum [a.u.]')

plt.tight_layout()
plt.show()
