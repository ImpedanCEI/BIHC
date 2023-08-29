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
from tqdm import tqdm

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

    def getShiftedPloss(self, Z, shift=20e6):
        '''
        Computes the power loss, shifting the impedance 
        curve rigidly in steps given by `shift`, to overlap
        with the spectral lines, giving a best (away from the
        line) and worst (on top of the line) case scenario.

        Parameters
        ----------
        Z : object
            Impedance object returned by Impedance class with 
            the frequency information of the impedance map given 
            by the user
        shift : float
            Frequency shift to be applied in steps to the 
            impedance curve
        '''

        [f,S] = self.spectrum
        deltaF = f[1]-f[0]
        fmax = Z.f[-1]

        #if impedance file is too short, we zero padd it
        #Otherwise the interpolation will assume for the 
        #missing frequencies a constant value 

        Zint = np.interp(f, Z.f, Z.Zr)
        Zmod = Z 
        Zmod.f, Zmod.Zr = f, Zint

        if Zmod.f[-1] > fmax: 
            mask = Zmod.f > fmax
            Zmod.Zr[mask] = 0.0

        size = int(shift/deltaF)
        shifts = np.arange(-size, size, 1, dtype=int) #shifting every frev

        power = np.array([])
        for step in tqdm(shifts, "Computing scan: ", total=2*size):
            Zmod.Zr = np.roll(Zint, step) 
            if step > 0: Zmod.Zr[:step] = 0.0
            if step < 0: Zmod.Zr[:2*size-step] = 0.0
            power = np.append(power, self.getPloss(Zmod)[0])

        return shifts, power

    def get2BeamPloss(self, Z_0, tau_s=None, s=None, offset1=None, offset2=None, Z_1=None):
        '''
        Computes the power loss for the two beams case
        given impedance object and the pahse_shift between 
        the two beams

        Implemented by Francesco Giordano

        Parameters
        ----------
        Z_0 : object
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

        e = 1.621e-19
        t = self.longitudinalProfile[0]
        f0 = 1/(t[-1]-t[0])
        [f,S] = self.spectrum
        Zreal_0 = Z_0.Zr
        Zf_0 = Z_0.f

        #Z_0
        if np.max(f)>np.max(Z_0.f):
            mask1=f>=0
            mask2=f<=np.max(Z_0.f)
            mask=mask1*mask2
            f=f[mask]
            Zreal_0=np.interp(f,Z_0.f,Z_0.Zr)
            S=S[mask]
            
        elif np.max(f)<np.max(Z_0.f):
            mask1=Z_0.f>=0
            mask2=Z_0.f<=np.max(f)
            mask=mask1*mask2

            Zf_0 = Z_0.f[mask]
            Zreal_0 = Z_0.Zr[mask]
            Zreal_1 = Z_1.Zr[mask]
            mask3 = f>=0
            f = f[mask3]
            S = S[mask3]
            
            Zreal_0 = np.interp(f,Zf_0,Zreal_0)

        #Z_1
        if Z_1 is not None:
            Zreal_1 = Z_1.Zr
            Zf_1 = Z_1.f

            if np.max(f)>np.max(Z_1.f):
                mask1=f>=0
                mask2=f<=np.max(Z_1.f)
                mask=mask1*mask2
                f=f[mask]
                Zreal_1=np.interp(f,Z_1.f,Z_1.Zr)
                S=S[mask]
                
            elif np.max(f)<np.max(Z_1.f):
                mask1=Z_0.f>=0
                mask2=Z_0.f<=np.max(f)
                mask=mask1*mask2

                Zf_1 = Z_1.f[mask]
                Zreal_1 = Z_1.Zr[mask]
                mask3 = f>=0
                f = f[mask3]
                S = S[mask3]

                Zreal_1 = np.interp(f,Zf_1,Zreal_1)
        
            mask3 = f>=0
            Zreal_1 = Zreal_1[mask3]
            Zf_1 = f

        mask3 = f>=0
        f = f[mask3]
        S = S[mask3]   
        Zreal_0 = Zreal_0[mask3]
        Zf_0 = f

        if Z_1 is not None and (offset1 is not None or offset2 is not None):

            if offset2 is None:
                offset2 = np.zeros_like(offset1)
            if offset1 is None:
                offset1 = np.zeros_like(offset2)

            # Formula with beam offsets
            P = (2*f0*e*self.filledSlots*self.Np*S)**2
            P_loss = []
            P_density_list = []
            for i, shift in tqdm(enumerate(tau_s), "Computing 2-beam power with offset: "):
                P_density = P*(Zreal_0+(offset1[i]+offset2[i])*Zreal_1)*(1-np.cos(2*np.pi*f*shift))
                P_loss.append(np.sum(P_density))  
        else:
            # Simplified formula
            P = (2*f0*e*self.filledSlots*self.Np*S)**2
            P_loss = []
            P_density_list = []
            for shift in tqdm(tau_s, "Computing 2-beam power: "):
                P_density = P*Zreal_0*(1-np.cos(2*np.pi*f*shift))
    #           P_density_list.append(P_density)
                P_loss.append(np.sum(P_density))  

        if self.verbose:
            pass
            #print(f'Computed Power loss: {P_loss} W') 

        self.P_loss = P_loss
        return P_loss
