'''
This example shows how to compute the dissipated power
for an SPS spectrum measurement

The impedance curve used is generated from a resonator

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


def fillingSchemeSPS_standard(ntrains):
    '''
    Returns the filling scheme for the SPS

    Parameters
    ----------
    ntrains: number of trains (batches)
    nbunches: default, 48. number of bunches per train
    '''
    # Define filling scheme: parameters
    nslots = 920 # Defining total number of slots for SPS
    nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72.. 
    batchspacing = 9 # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True]*nbunches
    st = [False]*batchspacing
    sc = [False]*(nslots - (nbunches+batchspacing)*ntrains)
    an = (bt + st)*ntrains + sc

    return an


# Set beam object with spectrum = 'user'

beam = bihc.Beam(Np=1.8e11, bunchLength=1.55e-9, fillingScheme=fillingSchemeSPS_standard(4), machine='SPS', spectrum='user') 

# Importing spectrum from measurement
import pandas as pd
def readFile(fileList): #--- read SA file - standard csv reader
    #df = pd.read_csv(path+fileList,  sep=',', engine='python')
    file=open(fileList, 'r')
    data=[row for row in csv.reader(file)]
    file.close()
    df0 = pd.DataFrame(data)
    index=df0.index[df0[0]=='DATA'][0]
    
    #for j in range(0,len(df0)):
    #    if df0[0][j] == 'DATA':
    #        index = j
    #        break
    hd=np.array(df0[:index])
    xd=np.array(df0[index+1:][0], dtype=float); 
    yd=np.array(df0[index+1:][1],dtype=float)
    return hd, xd, yd


# Set spectrum for beam
path = '/eos/user/a/avanel/Documents/MDshare/SPS_MD_04052023/Scan_2023-05-04_SA222/beam_spectra_LHC2_1e-05GHzTo2.0GHz_2023-05-04_11-13-06/'
header, xdata, ydata = readFile(path+'trace1.csv')
mask = xdata > 1

S = np.sqrt(10.0**((ydata[mask]-np.max(ydata[mask]))/10)) 
f = xdata[mask]
spectrum = [f, S]
beam.setSpectrum(spectrum) 

plt.plot(beam.spectrum[0], beam.spectrum[1])
plt.show()

# Obtain dissipated power from impedance resonator
Rs = 1.7e3
Qr = 250
fr = 0.786e9
Z = bihc.Impedance()
Z.getResonatorImpedance(Rs, Qr, fr)

# built-in plot spectrum and normalized impedance
beam.plotSpectrumAndImpedance(Z)

# Computing the dissipated power value
print(f'Beam 6675 power loss: {beam.getPloss(Z)[0]} W')
