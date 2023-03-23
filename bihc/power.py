'''
Power module to manage power loss computations
for one beam and two beams case.

Power loss is computed for a given beam (or beams)
power spectrum and the specified impedance map in
frequency.

* date: 12/12/2022
* author: Francesco Giordano, Elena de la Fuente, Leonardo Sito
'''

import numpy as np

class Power():
    '''Power Mixin class

    Class to encapsulate power computation methods
    It is inherited by Beam class

    Methods
    -------
    getPloss(Z)
        Computes the power loss for a given impedance object
    get2BeamPloss(Z, phase_shift)
        Computes the power loss for the 2 beams case for a 
        given impedance object and a phase shift between the
        beams
    '''
    def getPloss(self,Z):
        '''
        Computes the power loss for a given impedance object
        Implemented by Francesco Giordano

        Parameters
        ----------
        Z : object
            Impedance object returned by Impedance class with 
            the frequency information of the impedance map given 
            by the user
        '''

        e=1.621e-19
        t=self.longitudinalProfile[0]
        f0=1/(t[-1]-t[0])
        [f,S]=self.spectrum

        Zreal=Z.Zr
        Zf=Z.f

        if np.max(f)>np.max(Z.f):
            mask1=f>=0
            mask2=f<=np.max(Z.f)
            mask=mask1*mask2
            f=f[mask]
            Zreal=np.interp(f,Z.f,Z.Zr)
            S=S[mask]
            
        elif np.max(f)<np.max(Z.f):
            mask1=Z.f>=0
            mask2=Z.f<=np.max(f)
            mask=mask1*mask2
            Zf=Z.f[mask]
            Zreal=Z.Zr[mask]
            mask3= f>=0
            f=f[mask3]
            S=S[mask3]
            
            Zreal=np.interp(f,Zf,Zreal)
    
        mask3= f>=0
        f=f[mask3]
        S=S[mask3]   
        Zreal=Zreal[mask3]
        Zf=f

        P=f0*e*S*self.filledSlots*self.Np
        P_density=2*(P**2)*Zreal
        P_loss=np.sum(P_density) 
        if self.verbose:
            print(f'Computed Power loss: {P_loss} W') 

        return P_loss, P_density;

    def get2BeamPloss(self, Z, tau_s=None, s=None):
        '''
        Computes the power loss for the two beams case
        given impedance object and the pahse_shift between 
        the two beams

        Implemented by Francesco Giordano

        Parameters
        ----------
        Z : object
            Impedance object returned by Impedance class with 
            the frequency information of the impedance map given 
            by the user
        tau_s : float list
            Phase shift values between the two beams in seconds [s]
        s : float list
            Distances from the interaction point in [m]
        '''

        self.s = s
        if tau_s is not None:
            self.tau_s = tau_s
        elif s is not None:
            self.tau_s = 2*s/c
        else: 
            raise Exception('Specify s (distance from IP) or tau_s (phase shift of the two beams)')

        e=1.621e-19
        t=self.longitudinalProfile[0]
        f0=1/(t[-1]-t[0])
        [f,S]=self.spectrum
        Zreal=Z.Zr
        Zf=Z.f

        
        if np.max(f)>np.max(Z.f):
            mask1=f>=0
            mask2=f<=np.max(Z.f)
            mask=mask1*mask2
            f=f[mask]
            Zreal=np.interp(f,Z.f,Z.Zr)
            S=S[mask]
            
        elif np.max(f)<np.max(Z.f):
            mask1=Z.f>=0
            mask2=Z.f<=np.max(f)
            mask=mask1*mask2
            Zf=Z.f[mask]
            Zreal=Z.Zr[mask]
            mask3= f>=0
            f=f[mask3]
            S=S[mask3]
            
            Zreal=np.interp(f,Zf,Zreal)
    
        mask3= f>=0
        f=f[mask3]
        S=S[mask3]   
        Zreal=Zreal[mask3]
        Zf=f


        P=f0*e*S*self.filledSlots*self.Np
        P_loss = []
        P_density_list = []
        for shift in tau_s:
            P_density=4*(P**2)*Zreal*(1-np.cos(2*np.pi*f*shift))
#           P_density_list.append(P_density)
            P_loss.append(np.sum(P_density))  

        if self.verbose:
            print(f'Computed Power loss: {P_loss} W') 

        self.P_loss = P_loss
        return P_loss
