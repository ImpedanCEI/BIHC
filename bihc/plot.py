'''
Plot module to manage built-in plotting functions

It allows to easily plot filling schemes from timber,
or user defined beams. It has customized rcParams for
scientific plotting.

* date: 12/12/2022
* author: Francesco Giordano, Elena de la Fuente, Leonardo Sito
'''
import os, sys
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import rcParams, cycler

rcParams={
    # Set color cycle: blue, green, yellow, red, violet, gray
    'axes.prop_cycle' : cycler('color', ['0C5DA5', '00B945', 'FF9500', 'FF2C00', '845B97', '474747', '9e9e9e']),

    # Set default figure size
    'figure.figsize' : [4.55, 3.42],
    'figure.dpi': 160,
    'figure.autolayout': True,

    # Set x axis
    'xtick.direction' : 'in',
    'xtick.major.size' : 3,
    'xtick.major.width' : 0.5,
    'xtick.minor.size' : 1.5,
    'xtick.minor.width' : 0.5,
    'xtick.minor.visible' : True,
    'xtick.top' : True,

    # Set y axis
    'ytick.direction' : 'in',
    'ytick.major.size' : 3,
    'ytick.major.width' : 0.5,
    'ytick.minor.size' : 1.5,
    'ytick.minor.width' : 0.5,
    'ytick.minor.visible' : True,
    'ytick.right' : True,

    # Set line widths
    'axes.linewidth' : 0.5,
    'grid.linewidth' : 0.5,
    'lines.linewidth' : 1.,

    # Remove legend frame
    'legend.frameon' : False,

    # Always save as 'tight'
    'savefig.bbox' : 'tight',
    'savefig.pad_inches' : 0.05,

    # Use serif fonts
    # font.serif : Times,
    'font.family' : 'serif',
    'mathtext.fontset' : 'dejavuserif',
    }

def progressbar(it, prefix="", size=60, out=sys.stdout, count=None): # Python3.6+
    if count is None:
        count = len(it)
    def show(j):
        x = int(size*j/count)
        percent = int(j/count*100)
        #print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count}", end='\r', file=out, flush=True)
        print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {percent}%", end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)

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

        plt.title('Beam Longitudinal Profile')
        plt.ylim(0,np.max(s)/1e6)
        plt.xlim(tmin*1e6, tmax*1e6)            
        plt.tick_params(axis='both', which='major')
        plt.xlabel(r'time $[\mu s]$')
        plt.ylabel(r'Longitudinal time distribution $[1/\mu s]$')
        plt.legend(loc='upper right')
        plt.grid(True, color='gray', linestyle=':')
        plt.show()


    def plotPowerSpectrum(self, fmin=0, fmax=2, save=True, transparent=True):
        
        [f,S]=self.spectrum

        plt.figure()
        if (self._fillNumber!=0):
            label='fill: '+str(self._fillNumber)
        else:
            label=''

        plt.plot(f/1e9,S**2,label=label)

        plt.title('Beam Power Spectrum')
        plt.grid(True, color='gray', linestyle=':')
        plt.ylim(0,1)
        plt.xlim(fmin,fmax)
        plt.tick_params(axis='both', which='major')
        plt.xlabel("f [GHz]")
        plt.ylabel("Power Spectrum")

        if label:
            plt.legend()

        plt.show()

    def plotPowerSpectrumFromTimber(startDate,beamNumber):
        try:
            import pytimber
        except:
            print('This method uses pytimber. Please follow the installation guide to set it in your python environment')

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
        plt.grid(True, color='gray', linestyle=':')
        plt.tick_params(axis='both', which='major')
        plt.xlabel('f [GHz]')
        plt.ylabel('Power spectrum [a.u.]')
        plt.legend()
        if save: plt.savefig('PowerSpectrumFromTimber.png', transparent=True)
        plt.show()
        
        return f,S,timeStamps[-1]

    def plotImpedance(self,fMin=0,fMax=2e9):
        
        plt.figure()
        plt.plot(self.f/1e9,self.Zr,label='Re[Z]')
        plt.plot(self.f/1e9,self.Zi,color='r',label='Im[Z]',alpha=0.7)
        plt.tick_params(axis='both', which='major')
        plt.xlim(fMin/1e9,fMax/1e9)
        plt.ylim(0, )
        plt.title('Impedance')
        plt.xlabel('f [GHz]')
        plt.ylabel(r'Z $[\Omega]$')
        plt.grid(True, color='gray', linestyle=':')
        plt.legend()
        plt.show()

    def plotSpectrumAndImpedance(self, Z): #TODO: not normalize but have 2 y axis

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

        plt.figure()
        #spectrum
        plt.plot(f/1e9,S,label='Spectrum')
        ax = plt.gca()
        plt.xlim(0,self.fmax/1e9)
        plt.ylim(0,)
        plt.tick_params(axis='both', which='major')
        plt.xlabel("frequency [GHz]")
        plt.ylabel("Normalized Spectrum [a.u.]")
        plt.grid(True, color='gray', linestyle=':')
        plt.legend()
        #impedance
        axx = ax.twinx()
        axx.plot(f/1e9,Zreal, color='r', label='Impedance')
        axx.set_ylabel('Impedance [$\Omega$]')
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
            plt.grid(True, color='gray', linestyle=':')
            plt.ylim(0,)
            plt.xlim(np.min(t1)*1e6,np.max(t1)*1e6)
            plt.tick_params(axis='both', which='major')
            plt.xlabel("time [us]",fontsize=18)
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
        plt.grid(True, color='gray', linestyle=':')
        plt.ylim(0,)
        plt.xlim(np.min(t1)*1e6,np.max(t1)*1e6)
        plt.tick_params(axis='both', which='major')
        plt.xlabel("time [us]",fontsize=18)
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