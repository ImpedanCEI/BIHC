import numpy as np
import pandas as pd

from bihc.plot import Plot

class Impedance(Plot):
    
    def __init__(self, f, CST_file=0):
        self.isResonatorImpedance=False
        self.isRWImpedance=False
        self.f=f
        self.df = pd.DataFrame()
        self.df['f']=f
        self.df['Zr']=np.zeros(len(f))
        if CST_file!=0:
            self.getImpedanceFromCST(CST_file)
          
        
    def getResonatorImpedance(self,Rs,Qr,fr):
         
        f=self.f
        mask1= f==0
        f[mask1]=1e-5
        Z=Rs/(1+ 1j*Qr*(f/fr - fr/f))
        
        self.df=pd.DataFrame()
        self.df['Zr']=np.real(Z)
        self.df['Zi']=np.imag(Z)
        self.df['f']=f

        self.fr=fr
        self.Rs=Rs
        self.Qr=Qr
        self.isResonatorImpedance=True
        
        return Z.real;
        
    def getRWImpedance(self,L,b,sigma):
        Z0=376.73
        c=3e8
        f=self.f
        
        Z=L/(2*np.pi*b)*np.sqrt(Z0*np.abs(2*np.pi*f)/(2*c*sigma))
        
        self.real=np.real(Z)
        self.imag=np.imag(Z)

        self.b=b
        self.sigma=sigma
        self.L=L
        self.isRWImpedance=True
        
        return Z.real;
        
    def getImpedanceFromCST(self, path, unit='GHz', skip_header=2, skip_footer=0):
        
        data = np.genfromtxt(path, skip_header=skip_header, skip_footer=skip_footer)
        data = pd.DataFrame(data)
        # data.columns = ['f', 'Zr']
        data.columns = ['f', 'Zr', 'Zi']
        
        if unit == 'GHz':
            data['f']= data['f']*1e9
        elif unit == 'MHz':
            data['f']= data['f']*1e6
  
        self.df=data
    
    def getImpedanceFromPandas(self, path, unit='GHz'):
        
        data = pd.read_csv(path)
        
        data.columns = ['f', 'Zr']
        
        if unit == 'GHz':
            data['f']= data['f']*1e9
        elif unit == 'MHz':
            data['f']= data['f']*1e6
  
        self.df=data
        
        self.isResonatorImpedance=False
        
        return data
    