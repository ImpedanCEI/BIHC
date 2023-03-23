'''
Impedance module to manage impedance object creation from a txt file exported 
from CST database or thanks to (i) the resisitive wall impedance model (for a 
circular pipe) or (ii) broad-band resonator model.

The created impedance consists in a list of numpy.ndarray: the first element is
the frequency vector, the second is the impedance.

* date: 12/12/2022
* author: Francesco Giordano, Elena de la Fuente, Leonardo Sito
'''

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
    """   
    def __init__(self, f=np.linspace(0.1,2e9,int(1e5)), CST_file=None):
    
        self.f = f
        self.Zr = np.zeros(len(f))
        self.Zi = np.zeros(len(f))

        self.isResonatorImpedance = False
        self.isRWImpedance = False

        if CST_file is not None:
            self.getImpedanceFromCST(CST_file)          
        
    def getResonatorImpedance(self,Rs,Qr,fr,f=np.linspace(0.1,2e9,int(1e5))):
        """Creating the impedance curve from the broad-band resonator model.

        This methods creates an impedance curve using the broad-band resonator 
        model. It requires a shunt impedance value, a quality factor value, and
        the resonant frequency value. 

        Parameters
        ----------
        Rs : float
            Shunt impedance in the broad-band resonator model [Ohm].
        Qr : float
            Quality factor in the broad-band resonator model.
        fr : float
            Resonant frequency in the broad-band resonator model [Hz].
        
        Returns
        -------
        [f, Z] : numpy.ndarray list 
            Impedance curve. Returns the list of numpy arrays [frequency, 
            complex impedance].
        """
        f = self.f
        mask1 = f == 0
        f[mask1] = 1e-5
        Z = Rs/(1+ 1j*Qr*(f/fr - fr/f))
        
        self.Zr = np.real(Z)
        self.Zi = np.imag(Z)
        self.Z = Z
        self.f = f

        self.fr = fr
        self.Rs = Rs
        self.Qr = Qr
        self.isResonatorImpedance = True
        
        return [self.f, self.Zr+1j*self.Zi]
        
    def getRWImpedance(self, L ,b, sigma, f=np.linspace(0.1,2e9,int(1e5))):
        """Creating the impedance curve from the resistive wall impedance model.

        This methods creates an impedance curve using the resistive wall 
        impedance model for a resistive pipe of circular section with. 
        It considers the thick wall regime.

        Parameters
        ----------
        L : float
            Length of the pipe in the resistive wall impedance model [m].                
        b : float
            Pipe radius in the resisitve wall impedance model [m].
        sigma : float
            Electrical conductivity of the pipe in the resistive wall
            impedance model [Siemens/m].
        
        Returns
        -------
        [f, Z] : numpy.ndarray list 
            Impedance curve. Returns the list of numpy arrays [frequency, 
            complex impedance].
        """
        Z0 = 376.73
        c = 3e8
        f = self.f
        
        # !n.b. problem with the signum
        Z = (1+1j*np.sign(f))*L/(2*np.pi*b)*np.sqrt(Z0*np.abs(2*np.pi*f)/(2*c*sigma))
        
        self.Zr = np.real(Z)
        self.Zi = np.imag(Z)
        self.Z = Z

        self.b = b
        self.sigma = sigma
        self.L = L
        self.isRWImpedance = True
        
        return [self.f, self.Zr+1j*self.Zi]
        
    def getImpedanceFromCST(self, path, unit='GHz', skip_header=2, skip_footer=0):
        """Creating the impedance curve from CST file.

        This methods creates an impedance curve using a txt file, usually
        exported from CST.

        Parameters
        ----------
        path : str
            File name with extension.                
        unit : str, default "GHz"
            Units of the frquency array from the simulation. Either GHz or MHz.
        skip_header : int, default 2
            Line index to skip to at the beginning of the txt file.
        skip_footer : int, default 0
            Line index from which to skip to at the end of the txt file.
        
        Returns
        -------
        [f, Z] : numpy.ndarray list 
            Impedance curve. Returns the list of numpy arrays [frequency, 
            complex impedance].
        """   
        data = np.genfromtxt(path, skip_header=skip_header, skip_footer=skip_footer)
        
        if unit == 'GHz':
            self.f= data[:,0]*1e9
        elif unit == 'MHz':
            self.f= data[:,0]*1e6
  
        self.Zr = data[:,1]
        self.Zi = data[:,2]
        self.Z = self.Zr+1j*self.Zi

        return [self.f, self.Z]
    
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
    
