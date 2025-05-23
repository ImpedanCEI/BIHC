# copyright ################################# #
# This file is part of the BIHC Package.      #
# Copyright (c) CERN, 2024.                   #
# ########################################### #

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
import copy

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
    def __init__(self, f=np.linspace(0.1,2e9,int(1e5)), Z=None, CST_file=None):
    
        self.f = f
        self.Zr = np.zeros(len(f))
        self.Zi = np.zeros(len(f))

        if Z is not None:
            self.Zr = np.real(Z)
            self.Zi = np.imag(Z)
            if len(self.f) != len(self.Zr):
                print('[!] frequency array and impedance data have different lengths')

        self.isResonatorImpedance = False
        self.isRWImpedance = False

        if CST_file is not None:
            self.getImpedanceFromCST(CST_file)     

    def __add__(self, Zn):
        Z = Impedance(self.f)

        if not np.array_equal(self.f, Zn.f):
            Z.Zr = np.interp(self.f, Zn.f, Zn.Zr)
            Z.Zi = np.interp(self.f, Zn.f, Zn.Zi)
        else:
            Z.Zr = Zn.Zr
            Z.Zi = Zn.Zi

        Z.Zr += self.Zr
        Z.Zi += self.Zi   

        return Z
        
    def getResonatorImpedance(self,Rs,Qr,fr):
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

        self.fr = fr
        self.Rs = Rs
        self.Qr = Qr
        self.isResonatorImpedance = True

        f = self.f
        mask1 = f == 0
        f[mask1] = 1e-5
        Z = Rs/(1+ 1j*Qr*(f/fr - fr/f))
        
        self.Zr = np.real(Z)
        self.Zi = np.imag(Z)
        self.Z = Z
        
        return [self.f, self.Zr+1j*self.Zi]
        
    def getRWImpedance(self, L ,b, sigma):
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
        try:
            self.Zi = data[:,2]
        except:
            self.Zi = np.zeros_like(self.Zr)

        self.Z = self.Zr+1j*self.Zi

        return [self.f, self.Z]
    
    def getImpedanceFromFunc(self, func, f=None):
        '''Gets impedance from a python function
        in the form of Z = func(f)

        Ensures compatibility to Xwakes wake objects
        https://github.com/xsuite/xwakes/

        Example:
        -------
        >>> import xwakes as xw
        >>> wake = xw.WakeResonator(r=1e8, q=10, f_r=1e9,
                                    kind='longitudinal')
        >>> import bihc
        >>> Z = bihc.Impedance(f=np.linspace(0,1e9,10000))
        >>> Z.getImpedanceFromFunc(wake.components[0].impedance)
        '''

        if f is not None:
            self.f = f
        
        self.Z = func(f)
        self.Zr = np.real(self.Z)
        self.Zi = np.imag(self.Z)

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
        
        return data
    
    def copy(self):
        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj

    def getFrequencyRegions(self, vlines=None, figsize=[12,6], dpi=200):
        '''Interactively select points 
        of the impedance curve. The list of
        points will be returned and saved in 
        `self.freqregions`. How to use:
        - place cursor + Spacebar --> pick
        - Delete --> unpick
        - Enter --> End picking

        Returns
        -------
        freqregions: list
            np.array of frequencies of picked points [Hz]
        '''
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.plot(self.f*1e-9, self.Zr, color='k')

        if vlines is not None:
            ax.vlines(vlines*1e-9, ymin=np.min(self.Zr), ymax=np.max(self.Zr), color='k', linestyle='dashed', alpha=0.4)

        ax.set_ylabel('Impedance [$\Omega$]')
        ax.set_xlabel('Frequency [GHz]')
        ax.set_title('Spacebar --> pick | Delete --> unpick | Enter --> Finish', color='red')
        fig.tight_layout()

        #get frequencies 
        points = plt.ginput(n=0, timeout=0, mouse_add=None, mouse_pop=None, mouse_stop=None)
        self.freqregions = np.array(points)[:,0]*1e9
        plt.close()   

        return self.freqregions

