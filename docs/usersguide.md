# User's Manual

This section offers a description of the example suite available in BIHC package. This section asumes that the user has followed the Installation guide and has installed BIHC succesfully. 

<div class="warning">
<div class="admonition-title">
Warning
</div>
This documentation is currently under development
</div>

## User defined beam filling scheme

This example shows the definition of a filling  scheme as a list of True or Falses values. It is then used along some parameters of the bunches 
(Number of protons and slot space) to define a Beam object. Using some built in methods, the time distribution and beam spectrum are plotted. The manual plotting is also included below.
Finally, the power loss is calculated for a simple Pillbox impedance computed with `CST Stdio`[^1].

```python
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
custom_beam = bihc.Beam(bunchShape='GAUSSIAN', beamNumber=1, fillingScheme=an, Nb=Np, d=t0, verbose=False)

# built-in plotting
custom_beam.plotLongitudinalProfile()
custom_beam.plotPowerSpectrum()
```

The manual plotting can be done as well by retrieven the attributes of the beam class:

```python
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

```

## Timber filling scheme
This example shows how to compute the dissipated power
for a specific LHC fill Number. 

The filling scheme and the bunch parameters are obtained from a fill number through Timber database python interface `pytimber`.
The impedance curve used is from a generic Pillbox cavity simulated with CST studio and then exported in `.txt` format.

:::{admonition} Note
This example uses data from the CERN Timber database and requires access to NxCALS. 
See the [Installation guide](installation.md) for how to setup `pytimber` to acces the Timber database.
:::

```python
import bihc

# Data retrival from timber, with different bunch profile shapes
b_6675 = bihc.Beam(fillNumber=6675, bunchShape='GAUSSIAN')
# Exporting frequency array, needed for the impedance curve
[f, S] = b_6675.spectrum

# Importing an impedance curve
impedance_file = 'PillboxImpedance.txt'

Z = bihc.Impedance(f)
Z.getImpedanceFromCST(impedance_file)

# built-in plot spectrum and normalized impedance
b_6675.plotSpectrumAndImpedance(Z)

# Computing the dissipated power value
print(f'Beam 6675 power loss: {b_6675.getPloss(Z)[0]} W')

```
## LPC fillling schemes comparison

This example shows how to read filling schemes generated
by the LPC tool, in `.csv` format. It plots the comparison of three different
filling shcemes and computes the power loss for each of them, displaying it in the plot. 

:::{tip}
The LPC file is generated with the online `LPC tool` available in this [link](https://lpc.web.cern.ch/schemeEditor.html).
:::

```python
import matplotlib.pyplot as plt
import bihc

# LPC csv file names
# downloaded from LPC: https://lpc.web.cern.ch/cgi-bin/fillingSchemeTab.py
file1 = '25ns_2748b_2736_2258_2374_288bpi_12inj.csv'
file2 = '25ns_2374b_2361_1730_1773_236bpi_13inj_hybrid_2INDIV.csv' 
file3 = '8b4e_1972b_1967_1178_1886_224bpi_12inj.csv'

#------- Ploss calculation ----------

# Create beam object
bl = 1.2e-9         # bunch length [s]
Np = 2.3e11         # bunch intensity [protons/bunch]
bunchShape = 'GAUSSIAN'     # bunch profile shape in time 
fillMode = 'FLATTOP'        # Energy
fmax = 2e-9                 # Maximum frequency of the beam spectrum [Hz]

beamBcms1 = bihc.Beam(LPCfile=file1, Np=Np, bunchLength=bl, bunchShape=bunchShape, fillMode=fillMode, fmax=fmax)
beamBcms2 = bihc.Beam(LPCfile=file2,  Np=Np, bunchLength=bl, bunchShape=bunchShape, fillMode=fillMode,fmax=fmax)
beam8b4e = bihc.Beam(LPCfile=file3,  Np=Np, bunchLength=bl,  bunchShape=bunchShape, fillMode=fillMode, fmax=fmax)

# Importing an impedance curve
impedance_file = 'PillboxImpedance.txt'
Z = bihc.Impedance()
Z.getImpedanceFromCST(impedance_file)

# Computing the dissipated power value for the different filling schemes
print(f'25ns_2748b_2736_2258_2374_288bpi_12inj power loss: {beamBcms1.getPloss(Z)[0]} W')
print(f'25ns_2760b_2748_2494_2572_288bpi_13inj power loss: {beamBcms2.getPloss(Z)[0]} W')
print(f'8b4e_1972b_1967_1178_1886_224bpi_12inj power loss: {beam8b4e.getPloss(Z)[0]} W')

```

## Statistical min. max power loss
This section showcases the method `beam.getShiftedPloss()`. With this method, the impedance curve is shifted rigidly in frequency steps defined by the `shifts` parameter. This allows to consider uncertainties in the simulated impedance curve or in the beam spectral lines.
The maximum power will be obtained when the impedance peaks overlaps with one or more spectral lines. There is one example available for the SPS, using a user defined filling scheme (SPS standard 25ns, 4 batches). The other example is for the LHC, for which an LPC file is used as input.

### SPS example

Available in `examples/minMaxDissipatedPowerSPS.py`

```python
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
````

### LHC example

Available in `examples/minMaxDissipatedPower.py`

```python
import sys
sys.path.append('../../')

import bihc
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

#LPC csv file name
# downloaded from LPC: https://lpc.web.cern.ch/cgi-bin/fillingSchemeTab.py
file = '25ns_2760b_2748_2494_2572_288bpi_13inj.csv'

#------- Ploss calculation ----------

# Create beam object
bl = 1.2e-9                 # bunch length [s]
Np = 2.3e11                 # bunch intensity [protons/bunch]
bunchShape = 'GAUSSIAN'     # bunch profile shape in time 
fillMode = 'FLATTOP'        # Energy
fmax = 2e9                 # Maximum frequency of the beam spectrum [Hz]

beam = bihc.Beam(LPCfile=file, Np=Np, bunchLength=bl, bunchShape=bunchShape, 
                machine='LHC', fillMode=fillMode, spectrum='numeric', fmax=fmax)
[f,S] = beam.spectrum

# Create Impedance object
impedance_file = 'PillboxImpedance.txt'
Z = bihc.Impedance(f)
Z.getImpedanceFromCST(impedance_file)

# Get unshifted ploss 
ploss, ploss_density = beam.getPloss(Z) 

#---------------- Rigid Shift power loss ------------------------------
shift = 20e6  # distance between shift steps [Hz]
shifts, power = beam.getShiftedPloss(Z, shift=shift)

print(f'Minimum dissipated power: P_min = {np.min(power)}, at step {shifts[np.argmin(power)]}')
print(f'Maximum dissipated power: P_max = {np.max(power)}, at step {shifts[np.argmax(power)]}')
print(f'Average dissipated power: P_mean = {np.mean(power)}')
```

## Bunch profiles comparison

This example compares different bunch profile shapes using a same filling scheme given by an `LPC Tool`generated `.csv` file and the same impedance curve from a genereic Pillbox cavity computed with `CST Studio`.

It plots the impact of the different bunch shapes in the beam spectrum and computes the difference in power loss between the different bunch profiles used. 

:::{tip}
The LPC file is generated with the online `LPC tool` available in this [link](https://lpc.web.cern.ch/schemeEditor.html).
:::


```python
import matplotlib.pyplot as plt
import bihc

#LPC beam filling scheme 
file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'

# Data retrival from timber, with different bunch profile shapes
b_gauss= bihc.Beam(LPCfile=file, bunchShape='GAUSSIAN', verbose=False)
b_bin = bihc.Beam(LPCfile=file, bunchShape='BINOMIAL', verbose=False)
b_cos = bihc.Beam(LPCfile=file, bunchShape='COS2', verbose=False)

# Storing spectra 
[f_gauss, S_gauss] = b_gauss.spectrum
[f_bin, S_bin] = b_bin.spectrum
[f_cos, S_cos] = b_cos.spectrum

# Plotting 
fig, axs = plt.subplots(3,1)

axs[0].plot(f_gauss, S_gauss)
axs[1].plot(f_bin, S_bin)
axs[2].plot(f_cos, S_cos)

for ax in axs:
	ax.set_xlim(0, 2e9)
	ax.set_ylim(0, 1)
	ax.set_xlabel('frequency [Hz]')
	ax.set_ylabel('Spectrum [a.u.]')

# Importing an impedance curve
impedance_file = 'PillboxImpedance.txt'
Z = bihc.Impedance(f_gauss)
Z.getImpedanceFromCST(impedance_file)

# Computing the dissipated power value for the different Bunch Profiles
power_gauss = b_gauss.getPloss(Z)[0]
power_bin = b_bin.getPloss(Z)[0]
power_cos = b_cos.getPloss(Z)[0]

axs[0].text(0.3, 0.8, f'Gaussian Power loss: {round(power_gauss,2)} W', transform=axs[0].transAxes, color='tab:blue', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
axs[1].text(0.3, 0.8, f'Binomial Power loss: {round(power_bin,2)} W', transform=axs[1].transAxes, color='tab:blue', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))
axs[2].text(0.3, 0.8, f'Cosine Squared Power loss: {round(power_cos,2)} W', transform=axs[2].transAxes,color='tab:blue', weight='bold', bbox = dict(facecolor = 'white', alpha = 0.6))

print(f'Gaussian power loss: {power_gauss} W')
print(f'Binomial power loss: {power_bin} W')
print(f'Cosine Squared power loss: {power_cos} W')

plt.tight_layout()
plt.show()
```

## Analytical vs Numeric spectrum computation

This example compares the dissipated power difference obtained between 
using a numeric spactrum calculation or the analytic formula by C.Zannini
using a same filling scheme given by an `LPC Tool` generated `.csv` file and 
the same impedance curve from a simple Pillbox resonator computed with `CST Studio`.

It plots the impact of the spectrum type chosen: `numeric` or `analytic` in the 
beam spectrum and computes the difference in power loss.

```python
import matplotlib.pyplot as plt
import numpy as np
import bihc


# Beam data with different bunch profile shapes from LPC beam filling scheme 
file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
profiles = ['GAUSSIAN'] 
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
```


[^1]: CST Studio, [https://www.cst.com](https://www.cst.com)
