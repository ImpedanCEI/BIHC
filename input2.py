import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Beam new is the new version of Francesco's library. It has been modified.
from beam_new import Beam, Impedance

mpl.style.use(['ggplot'])

# Data retrival from timber
b_6675 = Beam(fillNumber=6675, bunchShape='BINOMIAL')
#b_6675=Beam( M=1320, bunchLength=4*0.05/3e8, d = 25e-9, phi=0, realMachineLength=True, Nb=0.6e11, bunchShape='GAUSSIAN', beamNumber=1, fillMode='FLATTOP', machine='LHC')
[f,S] = b_6675.spectrum

# Impedance plot
impedance_file = 'Impedance_file.txt'

Z = Impedance(f)
Z.getImpedanceFromCST(impedance_file)
Z.plot(fMax=Z.df.f.max())

# Spectrum plot 
[f,S] = b_6675.spectrum
plt.figure()
plt.plot(f/1e9,S, label='Spectrum')
plt.plot(Z.df.f/1e9, Z.df.Zr/Z.df.Zr.max(), label='Impedance Normalized')


# plt.plot(Zt.df.f/1e9, Zt.df.Zr/Zt.df.Zr.max(), label='FLATBOTTOM')
plt.xlabel('f [GHz]')
plt.xlim(0,2)
plt.ylabel('Spectrum [a.u.]')
plt.ylim(0,)
plt.legend()
plt.show()

# Ploss computation
print('Fill6675 BBLRC No Port  power loss: %.2f' %b_6675.getPloss(Z)[0], 'W')
