'''
In this example there is the definition of a filling scheme as a list of True or
Falses values.
It is then used along some parameters of the bunches (Number of protons and slot
space) to define a Beam object.
Finally, using some built in methods, the time distribution and beam spectrum
are plotted.

@date: Created on 08/12/2022
@author: lsito
'''

import bihc

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Define filling scheme: parameters
ninj = 10 # Defining number of injections
nslots = 3564 # Defining total number of slots for LHC
ntrain = 4 # Defining the number of trains
nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 

batchS = 7 # Batch spacing in 25 ns slots
injspacing = 37 # Injection spacing in 25 ns slots
BS = 200 # Batch spacing in ns

Np = 1.2e11 # Number of protons per bunch
t0 = 25e-9 # Slot space

# Defining the trains as lists of True/Falses
bt = [True]*nbunches
st = [False]*batchS
stt = [False]*injspacing
sc = [False]*(nslots-(ntrain*nbunches*ninj+((ntrain-1)*(batchS)*ninj)+((1)*injspacing*(ninj))))
an1 = bt+ st +bt+ st+ bt+ st+ bt+ stt
an = an1 * ninj + sc # This is the final true false sequence that is the beam distribution

# Data retrival from timber
custom_beam = bihc.Beam(bunchShape='GAUSSIAN', beamNumber=1, fillingScheme=an, Nb=Np, d=t0)

custom_beam.plotLongitudinalProfile()
custom_beam.plotPowerSpectrum()

