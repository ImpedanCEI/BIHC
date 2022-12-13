import numpy as np
import csv

from bihc.plot import Plot

class Impedance(Plot):
    """Define a Beam Coupling Impedance curve.

    It is possible to define an impedance curve using: (i) the broad-band 
    resonator model, (ii) the resistive wall impedance of a circular beam pipe, 
    and (iii) importing an external impedance text file (like the one that can
    be exported from CST).

    Mixin class for the Beam class.
    
    Parameters
    ----------
    f : numpy.ndarray, default np.linspace(0.1,2e9,int(1e5))
        Frequency array of the impedance curve.

        This is the numpy array containing the frequency point of the impedance
        curve. by default it is a linearly spaced vector going from 0.1 Hz to 2 
        GHz with 10'000 points.

    CST_file : str, deafult 0
        File name with extension of the txt file containing the impedance curve.
    
    Attributes
    ----------
    f : numpy.ndarray, default np.linspace(0.1,2e9,int(1e5))
        Frequency array of the impedance curve [Hz].
    Zr : numpy.ndarray
        Real part of the impedance in the speciefied frequency points [Ohm].
    Zi : numpy.ndarray
        Imaginary part of the impedance in the speciefied frequency points [Ohm].
    fr : float
        Resonant frequency in the broad-band resonator model [Hz].
    Rs : float
        Shunt impedance in the broad-band resonator model [Ohm].
    Qr : float
        Quality factor in the broad-band resonator model.
    b : float
        Pipe radius in the resisitve wall impedance model [m].
    sigma : float
        Electrical conductivity of the pipe in the resistive wall
        impedance model [Siemens/m].
    L : float
        Length of the pipe in the resistive wall impedance model [m].
    """
    
    def __init__(self, f=np.linspace(0.1,2e9,int(1e5)), CST_file=0):
    
        self.f=f
        self.Zr=np.zeros(len(f))
        self.Zi=np.zeros(len(f))

        self.isResonatorImpedance=False
        self.isRWImpedance=False

        if CST_file!=0:
            self.getImpedanceFromCST(CST_file)          
        
    def getResonatorImpedance(self,Rs,Qr,fr):
         
        f=self.f
        mask1= f==0
        f[mask1]=1e-5
        Z=Rs/(1+ 1j*Qr*(f/fr - fr/f))
        
        self.Zr=np.real(Z)
        self.Zi=np.imag(Z)
        self.f=f

        self.fr=fr
        self.Rs=Rs
        self.Qr=Qr
        self.isResonatorImpedance=True
        
        return Z.real;
        
    def getRWImpedance(self,L,b,sigma):
        Z0=376.73
        c=3e8
        f=self.f
        
        # !n.b. problem with the signum
        Z=(1+1j*np.sign(f))*L/(2*np.pi*b)*np.sqrt(Z0*np.abs(2*np.pi*f)/(2*c*sigma))
        
        self.Zr=np.real(Z)
        self.Zi=np.imag(Z)

        self.b=b
        self.sigma=sigma
        self.L=L
        self.isRWImpedance=True
        
        return Z.real;
        
    def getImpedanceFromCST(self, path, unit='GHz', skip_header=2, skip_footer=0):
        
        data = np.genfromtxt(path, skip_header=skip_header, skip_footer=skip_footer)
        
        if unit == 'GHz':
            self.f= data[:,0]*1e9
        elif unit == 'MHz':
            self.f= data[:,0]*1e6
  
        self.Zr = data[:,1]
        self.Zi = data[:,2]
    
    def getImpedanceFromPandas(self, path, unit='GHz'):
        import pandas as pd

        data = pd.read_csv(path)
        data.columns = ['f', 'Zr']
        
        if unit == 'GHz':
            data['f']= data['f']*1e9
        elif unit == 'MHz':
            data['f']= data['f']*1e6
  
        self.f = data['f']
        self.Zr = data['Zr']
        self.isResonatorImpedance=False
        
        return data
    
