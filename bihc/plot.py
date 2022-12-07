import matplotlib.pyplot as plt 
import numpy as np
import pytimber

#customize plots
from matplotlib import rcParams, cycler
plt.style.use('ggplot')
rcParams['lines.linewidth'] = 2
rcParams['xtick.labelsize'] = 'xx-large'
rcParams['ytick.labelsize'] = 'xx-large'

class Plot():

    def plotLongitudinalProfile(self, tmin=-1, tmax=-1):

        [t,s]=self.longitudinalProfile

        if tmin==-1:
            tmin=np.min(t)
        if tmax==-1:
            tmax=np.max(t)

        label=str(self.filledSlots) + ' bunches'   
        
        plt.figure()
        plt.plot(t*1e6,s/1e6,label=label)
        plt.grid(True)
        plt.ylim(0,np.max(s)/1e6)
        plt.xlim(tmin*1e6, tmax*1e6)
            
        plt.tick_params(axis='both', which='major')
        plt.xlabel(r't $[\mu s]$')
        plt.ylabel(r'Longitudinal time distribution $[1/\mu s]$')
        plt.legend(loc='upper right')
        plt.show()


    def plotPowerSpectrum(self, fmin=0, fmax=2, save=True, transparent=True):

        if (self._fillNumber!=0):
            label='fill: '+str(self._fillNumber)
        else:
            label=''
        
        [f,S]=self.spectrum
        
        plt.figure()
        plt.plot(f/1e9,S**2,label=label)
        plt.grid('on')
        plt.ylim(0,)
        plt.xlim(fmin,fmax)
        plt.tick_params(axis='both', which='major')
        plt.xlabel("f [GHz]")
        plt.ylabel("Power Spectrum")
        plt.legend()
        plt.show()


    def plotSpectrumFromTimber(startDate,beamNumber):
        db=pytimber.LoggingDB()
        
        if(beamNumber==1):
            freq='ALB.SR4.B1:SPEC_FREQ'
            powSpec='ALB.SR4.B1:SPEC_POW'
        elif(beamNumber==2):
            freq='ALB.SR4.B2:SPEC_FREQ'
            powSpec='ALB.SR4.B2:SPEC_POW'  
    
        delayTime=5*60
        
        t1 = pytimber.parsedate(startDate) - delayTime
        t2= 'next'
        
        f=db.get(freq,t1,t2)  #Siccome i filledbuckets non cambiano prendo tutti i valori dall'inizio 
                                         #al FLATTOP e poi seleziono solo l'ultimo vettore
        timeStamps,f=f[freq]
        f=f[-1]/1e9
    
        S=db.get(powSpec,t1,t2)
        timeStamps,S=S[powSpec]
        print('PowerSpectrum Date :', pytimber.dumpdate(timeStamps[-1]))
        S=S[-1]
        plt.figure()
        plt.plot(f, S, color='b',label=pytimber.dumpdate(timeStamps[-1]))
        ply.xlim(0,2)
        plt.ylim(np.min(S),np.max(S))
        plt.grid('on')
        plt.tick_params(axis='both', which='major')
        plt.xlabel('f [GHz]')
        plt.ylabel('Power spectrum [a.u.]')
        plt.legend()
        if save: plt.savefig('PowerSpectrumFromTimber.png', transparent=True)
        plt.show()
        
        return f,S,timeStamps[-1]

    def plotImpedance(self,fMin=0,fMax=2e9):
        
        plt.figure()
        plt.plot(self.df['f']/1e9,self.df['Zr'],label='Re[Z]')
        plt.plot(self.df['f']/1e9,self.df['Zi'],color='r',label='Im[Z]',alpha=0.7)
        plt.tick_params(axis='both', which='major')
        plt.xlim(fMin/1e9,fMax/1e9)
        plt.ylim(0, )
        plt.title('Impedance')
        plt.xlabel('f [GHz]')
        plt.ylabel(r'Z $[\Omega]$')
        plt.grid(True)
        plt.legend()
        plt.show()

    def plotSpectrumAndImpedance(self, Z):

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

        plt.figure()
        plt.plot(f/1e9,S,label='normalized spectrum')
        plt.plot(f/1e9,Zreal, color='r', label='Impedance [ohm]')
        plt.xlim(0,2)
        plt.ylim(0,)
        plt.tick_params(axis='both', which='major')
        plt.xlabel("f [GHz]")
        plt.ylabel("Impedance")
        plt.grid('on')
        plt.legend()
        plt.show()

    def plotCollide(beam1, beam2):
        
        [t1,s1]=beam1.longitudinalProfile
        [t2,s2]=beam2.longitudinalProfile
        s2=s2[::-1]
        step=20000
                    
        IP=(np.max(t1)/2)
        deltaT=(np.max(t1)/7)
        mask1=t1<IP + deltaT
        mask2=t1>IP - deltaT
        mask=mask1*mask2
        for j in range(20):
            s1=np.roll(s1,step)
            s2=np.roll(s2,-step)
            
            plt.figure()
            plt.plot(t1[mask]*1e6,s1[mask],label='beam1',color='b')
            plt.plot(t2[mask]*1e6,s2[mask],label='beam2',color='r',alpha=0.7)
           #plt.plot(t2[mask]*1e6,(s2+s1)[mask],label='beam1+beam2',color='g',alpha=0.7)
            plt.grid('on')
            plt.ylim(0,)
            plt.xlim(np.min(t1)*1e6,np.max(t1)*1e6)
            plt.tick_params(axis='both', which='major')
            plt.xlabel("t [us]",fontsize=18)
            plt.ylabel("Longitudinal Profile [a.u.]")
            plt.legend()
            plt.show()

    def plot2Beam(beam1, beam2, shift=0):

        [t1,s1]=beam1.longitudinalProfile
        [t2,s2]=beam2.longitudinalProfile
        deltaT=t1[1]-t1[0]
        step=int(shift/deltaT)
        s1=np.roll(s1,step)
        s2=np.roll(s2,-step)

        plt.figure()
        plt.plot(t1*1e6,s1,label='beam1',color='b')
        plt.plot(t2*1e6,s2,label='beam2',color='r',alpha=0.7)
#        plt.plot(t2*1e6,(s2+s1)[mask],label='beam1+beam2',color='g',alpha=0.7)
        plt.grid('on')
        plt.ylim(0,)
        plt.xlim(np.min(t1)*1e6,np.max(t1)*1e6)
        plt.tick_params(axis='both', which='major')
        plt.xlabel("t [us]",fontsize=18)
        plt.ylabel("Longitudinal Profile [a.u.]")
        plt.legend(fontsize=16)
        plt.show()

    def saveLongitudinalDistribution(self, path, phase_shift=0, normalise=False, plot=False):

        #phase_shift is intendet to be in ns 
        [t,s]=self.longitudinalProfile
        phase_shift = phase_shift*1e-9
        dt = t[1] - t[0]
        roll_points = phase_shift/dt
        s=np.roll(s,int(roll_points))
        if normalise:
            s = s/np.max(s)
        data = np.array([t, s])
        data = data.T
        np.savetxt(fname=path ,X=data, header='t[s]:     s[a.u]:',fmt=['%e','%e'], delimiter='     ')