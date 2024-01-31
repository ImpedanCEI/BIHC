'''
# Power loss computations with BIHC

- code: https://github.com/lsito/BIHC

- docs: https://bihc.readthedocs.io/en/latest/

- install: `pip install bihc`

- upgrade: `pip install bihc --upgrade`
'''

import numpy as np
import matplotlib.pyplot as plt
import bihc

# LPC generated filling scheme file name
# ---------------------------------
LPCfile ='25ns_2760b_2748_2494_2572_288bpi_13inj.csv' #fill from LPCTool website
Np = 1.63e11    # Intensity p/b
bl = 1.2e-9     # bunch length
fmax = 2.6790152140831e9 #from impedance file
      
# Defining beam object
# ---------------------------------
beam = bihc.Beam(Np=Np, bunchLength=bl, LPCfile=LPCfile, bunchShape='GAUSSIAN', fmax=fmax) 

# Defining impedance
# ---------------------------------
file = 'impedances/63wcase.txt' 
Z = bihc.Impedance()
Z.getImpedanceFromCST(file)
Z.Zr = np.abs(Z.Zr)

# Compute power loss
# ---------------------------------
power1b = beam.getPloss(Z)[0]
shifts, power1b_s = beam.getShiftedPloss(Z, shift=20e6) 
#gives worst case by shifting the impedance curve
#until on top 

# Plot impedance and power spectrum
# ---------------------------------
fig, ax = plt.subplots()

axx = ax.twinx()

l1, = ax.plot(beam.powerSpectrum[0], beam.powerSpectrum[1], color='b', alpha=0.8)
l2, = axx.plot(Z.f, Z.Zr, color='k')

# Shifted impedance curve
step = shifts[np.argmax(power1b_s)]
Zshifted = np.roll(Z.Zr, step)
if step > 0: Zshifted[:step] = Z.Zr[0]
l3, = axx.plot(Z.f, Zshifted, color='r', alpha=0.8)

ax.legend([l1, l2, l3], ['$\Lambda^2$', file.split('/')[1].split('.txt')[0], file.split('/')[1].split('.txt')[0]+' shifted'], loc=2)
ax.set_ylabel('Power spectrum amplitude')
axx.set_ylabel('Real Impedance [Ohm]')
ax.set_xlabel('Frequency [Hz]')
ax.set_xlim((0, fmax))
ax.set_ylim(ymin=0)
axx.set_ylim(ymin=0)
ax.set_title('Impedance and power spectrum')
ax.text(0.95, 0.1, f'Estimated power: {round(power1b,2)} W', ha='right', va='center', transform=ax.transAxes)
ax.text(0.95, 0.05, f'Max. Estimated power: {round(np.max(power1b_s),2)} W', c='r',ha='right', va='center', transform=ax.transAxes)
plt.show()

