import numpy as np

class Power():

    def getPloss(self,Z):

        if(self._NbIsComputed==False):
            self.setNbFromFillNumber()

        e=1.621e-19
        t=self.longitudinalProfile[0]
        f0=1/(t[-1]-t[0])
        [f,S]=self.spectrum

        Zreal=Z.df['Zr']
        Zf=Z.df['f']

        if np.max(f)>np.max(Z.df['f']):
            mask1=f>=0
            mask2=f<=np.max(Z.df['f'])
            mask=mask1*mask2
            f=f[mask]
            Zreal=np.interp(f,Z.df['f'],Z.df['Zr'])
            S=S[mask]
            
        elif np.max(f)<np.max(Z.df['f']):
            mask1=Z.df['f']>=0
            mask2=Z.df['f']<=np.max(f)
            mask=mask1*mask2
            Zf=Z.df['f'][mask]
            Zreal=Z.df['Zr'][mask]
            mask3= f>=0
            f=f[mask3]
            S=S[mask3]
            
            Zreal=np.interp(f,Zf,Zreal)
    
        mask3= f>=0
        f=f[mask3]
        S=S[mask3]   
        Zreal=Zreal[mask3]
        Zf=f

        P=f0*e*S*self.filledSlots*self.Nb
        P_density=2*(P**2)*Zreal
        P_loss=np.sum(P_density)  
        return P_loss, P_density;

    def get2BeamPloss(self, Z, phase_shift):

        e=1.621e-19
        t=self.longitudinalProfile[0]
        f0=1/(t[-1]-t[0])
        [f,S]=self.spectrum
        Zreal=Z.df['Zr']
        Zf=Z.df['f']
        
        if np.max(f)>np.max(Z.df['f']):
            mask1=f>=0
            mask2=f<=np.max(Z.df['f'])
            mask=mask1*mask2
            f=f[mask]
            Zreal=np.interp(f,Z.df['f'],Z.df['Zr'])
            S=S[mask]
            
        elif np.max(f)<np.max(Z.df['f']):
            mask1=Z.df['f']>=0
            mask2=Z.df['f']<=np.max(f)
            mask=mask1*mask2
            Zf=Z.df['f'][mask]
            Zreal=Z.df['Zr'][mask]
            mask3= f>=0
            f=f[mask3]
            S=S[mask3]
            
            Zreal=np.interp(f,Zf,Zreal)
    
        mask3= f>=0
        f=f[mask3]
        S=S[mask3]   
        Zreal=Zreal[mask3]
        Zf=f


        P=f0*e*S*self.filledSlots*self.Nb
        P_loss = []
        P_density_list = []
        for shift in phase_shift:
            P_density=4*(P**2)*Zreal*(1-np.cos(2*np.pi*f*shift))
#           P_density_list.append(P_density)
            P_loss.append(np.sum(P_density))  

        return P_loss;# f, P_density_list;