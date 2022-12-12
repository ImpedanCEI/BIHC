'''
This example shows how to produce some impedance curves using some built in
methods. This is an alternative to importing an impedance curve from CST.


@date: Created on 09/12/2022
@author: lsito
'''

#TODO Adding more resontors impedances and the impedance from 
import bihc

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Data retrival from timber, with different bunch profile shapes
Z_RW = bihc.Impedance()
Z_R = bihc.Impedance()

Z_RW.getRWImpedance(L=1, b=12e-3, sigma=2.5e7)
Z_R.getResonatorImpedance(Rs=1, Qr=50, fr=1e9)

# Automatic plotting
Z_RW.plotImpedance()
Z_R.plotImpedance()

# Manual plotting
plt.figure()
plt.plot(Z_RW.f, Z_RW.Zr)
plt.plot(Z_RW.f, Z_RW.Zi)
plt.plot(Z_R.f, Z_R.Zr)
plt.plot(Z_R.f, Z_R.Zi)
plt.show()
