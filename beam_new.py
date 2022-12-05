import numpy as np
import pytimber
import pylab as pl
import sys
import pandas as pd

# Fancier way to show graphs
import matplotlib.pyplot as plt
plt.style.use('ggplot')

#%%
# There should be no need for (object), it is from Python 2
class Beam(object):
    '''Defining the beam characteristics.

    Enables to create a beam with certain characteristics. The beam can be
    created from an existing fill stored in Timber or from a self defined
    filling scheme. It is possible to assign the bunch lenght, the bunch
    shape, the total intensity.
    '''

    def __init__(self, M=3564, A=1, fillNumber=0,
                 bunchLength=1.2e-9, phi=0, realMachineLength=True,
                 ppbk=250,d=0 , Nb=2.3e11, bunchShape='GAUSSIAN', beamNumber=1, 
                 fillMode='FLATTOP', fillingScheme=[False]*3564, machine='LHC'):
        
        c = 299792458 # Speed of light in vacuum [m/s]
        
        self.M = M # Default max numebr of buckets
        self.A_GLOBAL = A
        self.BUNCH_LENGTH_GLOBAL = bunchLength/4 # Bunch lenght [s]
        
        self.PHI_GLOBAL = phi
        self.realMachineLength = realMachineLength
        self.ppbk = ppbk

        self._bunchShape = bunchShape # Bunche shape (analytical function)
        self.J = 1
        self.Nb = Nb # Number of particles per bunch
        self.fillMode = fillMode
        self._fillNumber = fillNumber # Fill number relative to a particular fill of the machine
        
         
        
        self._isSpectrumReady = False
        self.isATimberFill    = False
         
        self._machine=machine # Select the machine you are working with
        # Some parameters are set in a different way depending on the machine you are working with
        if self._machine== 'LHC':
            self.BUCKET_MAX=3564
            RING_CIRCUMFERENCE = 26658.883                 #[m]
            GAMMA_R = 7461
            
            
        elif self._machine=='SPS':
            pass
        elif self._machine== 'PS':
            self.BUCKET_MAX=21
            RING_CIRCUMFERENCE=628
            GAMMA_R = 27.7366 #28.7185  # p=26GeV
            
            
        if self.M > self.BUCKET_MAX:
            self.M = self.BUCKET_MAX
            print('Number of bucket that could be filled set to: ', self.M)
        
        BETA_R = np.sqrt(1 - (1/GAMMA_R**2))
        self.T_1_TURN = RING_CIRCUMFERENCE/(c*BETA_R)

        # d is the time (space) of one bucket
        if d==0:
            self.d = self.T_1_TURN/self.BUCKET_MAX
        else:
            self.d = d
        
        self.l = self.d/2
        self._bunchLength = np.zeros(self.M)
        self.phi = np.zeros(self.M)
        self._spectrum = [np.zeros(self.M), np.zeros(self.M)]    #[f,S]
        self._fillingScheme=fillingScheme[0:self.M]
        
        
        
        self._beamNumber = beamNumber # Beam number, either 1 or 2
        
        # Se hai dato un fillNumber, cioé è diverso da 0 che è il default allora chiama
        # il metodo setBeamFromFillNumber, altrimenti, se avevi dato un filling scheme,
        # e questo significa che il vettore non è di tutti 0, e controllo sommando tutti 
        # i valori tra di loro, allora chiama setCustomBeamWithFillingScheme, infine, terza
        # opzione è il setCustomBeam
        if fillNumber > 0:
            self.setBeamFromFillNumber(fillNumber, fillMode, beamNumber)
        elif sum(self._fillingScheme) > 0:
            self.setCustomBeamWithFillingScheme()
        else:
            self.setCustomBeam()
        
#%%     
    @property
    def bunchLength(self):
        return self._bunchLength
    
    @bunchLength.setter
    def bunchLength(self,newBunchLength):
        print('updating bunch length')
        self._bunchLength=np.zeros(self.M)
        self._bunchLength[self._fillingScheme]=newBunchLength/4
        self._setBunches()

#%%        
    @property
    def bunchShape(self):
        return self._bunchShape
    
    @bunchShape.setter
    def bunchShape(self,newShape):
        shapeList=['GAUSSIAN', 'BINOMIAL' , 'COS2']
        print('updating bunch shape...')
        if newShape in shapeList:
            self._bunchShape=newShape
            self._setBunches()
            self._isSpectrumReady=False
        else:
            raise ValueError("bunchShape should be: 'GAUSSIAN', 'BINOMIAL' or 'COS2'." )

#%%     
    
    
    @property
    def fillNumber(self):
        return self._fillNumber
    
    @fillNumber.setter
    def fillNumber(self, newFillNumber):
        print('updating fill number...')
        if newFillNumber > 0:
            self._fillNumber=newFillNumber
            self.setBeamFromFillNumber(self._fillNumber)
        else:
            raise ValueError("fillNumber should be a positive value")

#%%            
    
    @property
    def spectrum(self):
        if self._isSpectrumReady:
            return self._spectrum
        else:
            [t,s]=self.longitudinalProfile
            deltaT=t[10]-t[9] 
            fc=1/deltaT                                        #in frequency we have a periodic signal of period fc where fc is 1/deltaT where the sampling step
                                                               #We are interested only in the range between -fc/2 and fc/2 ( In particular (0,fc/2) beacouse x is real)
    
            S=np.fft.fft(s,len(s))
            S=np.fft.fftshift(S)
            S=S*deltaT
            deltaF=fc/len(s)
    
            print ('DC component: ', np.max(np.abs(S)))  
            f = np.linspace(-fc/2, fc/2 - deltaF, len(S))      #vector of K point from min_value to max_value (Domain in frequency)
            self._spectrum=[f, np.abs(S)]
            self._isSpectrumReady=True
            
            return self._spectrum
        
    @spectrum.setter
    def spectrum(self, newSpectrum):
#        self._spectrum=newSpectrum
        raise Exception("Spectrum can not be assigned")
                       

        
#%%
    
    def _setBunches(self):

        if((self.realMachineLength) and (self.d*self.M>self.T_1_TURN)):
            self.d=self.T_1_TURN/self.M
            print("d has been resized to ", self.d ," because it was to big to fill ",self.M," buckets in the real machine length")
        
        deltaD=self.d/self.ppbk
        tTemp = np.arange(-self.d/2, self.d/2 ,deltaD)
        
        print ("Elaborating Data...")
        
        s = np.zeros(len(tTemp)) #[0] * len(t)
        sTemp=np.zeros(len(s))
        H=0
        self.filledSlots=0
        
        for i in range(self.M):
            if self._fillingScheme[i]==True:
                
                self.filledSlots+=1
                if(self._bunchShape=='BINOMIAL'):
                    H=(np.sqrt(2/np.log(2)))*self._bunchLength[i]*4
                    
                    sTemp=( 1- 4*((tTemp - self.phi[i])/(H))**2) #Binomial function
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(sTemp>0)
                    sTemp=sTemp*mask2*mask
                    sTemp=2*(sTemp**2.5)/H
                    
                elif(self._bunchShape=='GAUSSIAN'):
                    sTemp=1/(self._bunchLength[i] * np.sqrt(2 * np.pi)) *(np.e)**(-((tTemp - self.phi[i])**2)/(2*self._bunchLength[i]**2))  #Gaussian function
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=np.ones(len(sTemp))
                    sTemp=sTemp*mask2*mask
                    
                elif(self._bunchShape=='COS2'):
                    tc=self._bunchLength[i]*2.77
                    sTemp=1/tc*(np.cos(np.pi*(tTemp - self.phi[i])/(2*tc)))**2  #cos^2 function
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(abs(tTemp)<tc)
                    sTemp=sTemp*mask2*mask
                elif(self._bunchShape=='PARABOLIC'):
                    
                    sTemp=(1-(1/(self._bunchLength[i]**2))*(tTemp- self.phi[i])**2)  #Parabolic function 
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(sTemp>0)       
                    sTemp=sTemp*mask*mask2
                    sTemp=sTemp/(sum(sTemp)*deltaD)
                elif(self._bunchShape=='q-GAUSSIAN'):
                    
                    H = 2 * self._bunchLength[i] * np.sqrt(2*np.log(2)/(1 - 2**(-0.4))) 
                    
                    sTemp=( 1- 4*((tTemp - self.phi[i])/(H))**2) #Binomial function
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(sTemp>0)
                    sTemp=sTemp*mask2*mask
                    sTemp=32*(sTemp**2.5)/(5 * np.pi * H)
               
                   #set to zero the part of each bunch that is outside the l range arount his mean
                
            
            else:
                sTemp=np.zeros(len(tTemp))
            if(i==0):
                s=sTemp
            else:
                s=np.concatenate((s,sTemp),axis=0)
        #ppbk sono i punti di un bucket, M*ppbk restituisce il numero di punti necessari per M buckets
        
        
        
        t_max=self.d*self.M-self.d/2
        
        if(self.realMachineLength):
            t=np.arange(-self.d/2,self.T_1_TURN -self.d/2,deltaD)     
        else:
            t=np.arange(-self.d/2,t_max,deltaD)
         
        s_end=np.zeros(len(t) - len(s))
        s=np.concatenate((s,s_end),axis=0)
        sTemp=s
        for k in range(1,self.J):
            s=np.concatenate((s,sTemp),axis=0)
        
        tn=t
        t=np.linspace(np.min(tn),np.max(tn),len(s))
        
        
        
        self.totalBeamCharge=self.Nb* 1.602176e-19*self.filledSlots

        s=s/(self.filledSlots*self.J)
        
        dt=t[1]-t[0]
    
        s_integr=np.sum(s)*dt
    
        mask=self._bunchLength!=0
        print ('Average bunch length = ', np.mean(self._bunchLength[mask])*4, 's')
        print ('Integral of s = ' , s_integr)
        print ("Buckets filled : ", self.filledSlots)
        print ('Nb: ',self.Nb)
        print ("Total beam charge: ", self.totalBeamCharge, "C")
        
        self.longitudinalProfile=[t,s]
        
#%%        

    def setBeamFromFillNumber(self,fillNumber,fillMode='FLATTOP',beamNumber=1):
        self._fillNumber=fillNumber
        self.fillMode=fillMode
        self._beamNumber=beamNumber
        self.isATimberFill=True

        print ('Downloading data from Timber...')
        db=pytimber.LoggingDB()        
        bunchLengths='LHC.BQM.B'+str(beamNumber)+':BUNCH_LENGTHS'
        filledBuckets='LHC.BQM.B'+str(beamNumber)+':FILLED_BUCKETS'

    
        fill=db.getLHCFillData(fillNumber)
        ts = fill['startTime']
#        print 'MODE of the selected fill: '
        for j in range(len(fill['beamModes'])):
#            print fill['beamModes'][j]['mode']
            if(fill['beamModes'][j]['mode']==fillMode):
                MODE=j
            
        t1=fill['beamModes'][MODE]['startTime']
        t2=fill['beamModes'][MODE]['endTime']
        print ('Mode selected: ',self.fillMode, 'starts at: ', pytimber.dumpdate(t1), 'ends at',pytimber.dumpdate(t2))
    
        fb=db.get(filledBuckets,ts,t2)
        timeStamps,fb=fb[filledBuckets]
  
        std=db.get(bunchLengths,t1,t2)
        timeStamps,std=std[bunchLengths]
        '''il ciclo while selezione l'ultimo vettore std e fb che non e' zero tra tutti quelli letti'''
        
        j=len(fb)-1
        while True:
            if np.sum(fb[j])!=0:
                break
            j=j-1
        fb=fb[j]
    
        i=len(std)-1
        while True:
            if np.sum(std[i])!=0:
                break
            i=i-1
        print ('Date of the loaded data:', pytimber.dumpdate(timeStamps[i]))
        std=std[i]
        
        self.beamDate=pytimber.dumpdate(timeStamps[i])
        self._bunchLength=np.zeros(len(std))
        self.phi=np.zeros(len(std))
        FB=np.zeros(len(fb))
    
    
        for j in range(len(fb)):
            good=isinstance( fb[j], (float) )
            if(good and fb[j]!=0):
                FB[j]=int((fb[j] -1)/10)
            else:
                FB[j]=-1
            
        for j in range(len(std)):    #std is a sorted vector where std[i[j]] rappresent the correct position in the time flow
            stdIsGood=isinstance( std[j], (float) )
            if((stdIsGood) and (FB[j]!=-1) and (FB[j]<self.M)):
                self._bunchLength[int(FB[j])]=std[j]/4
                self._fillingScheme[int(FB[j])]=True
         

        self._bunchLength=self._bunchLength[0:self.M]
        self.phi=self.phi[0:self.M]
        self.setNbFromFillNumber()
        self._setBunches()
        

#%%        
        
    def setCustomBeam(self):

        self.isATimberFill=False

        self._fillingScheme[0:self.M]=[True]*self.M
        self._bunchLength[self._fillingScheme]=self.BUNCH_LENGTH_GLOBAL  #std vector of a single turn in the machine
#        self._fillingScheme[:M-1]=True
#        self.A[self._fillingScheme]=self.A_GLOBAL   #amplitude vector of a single turn in the machine
        self.phi=np.ones(self.M)*self.PHI_GLOBAL
            
        self._setBunches()

#%%        
        
    def setCustomBeamWithFillingScheme(self):
        print('Setting custom beam from filling scheme')
        if(len(self._fillingScheme)>self.M):
            sys.exit("ERROR : the length of the fillingScheme exceed the slot that you set (M)")
        elif(len(self._fillingScheme)<self.M):
            padding=[False]*(self.M - len(self._fillingScheme))
            self._fillingScheme=np.concatenate((self._fillingScheme,padding),axis=0)

        
        self._bunchLength=np.zeros(self.M)  #std vector of a single turn in the machine
        self.A=np.zeros(self.M)   #amplitude vector of a single turn in the machine
        self.phi=np.zeros(self.M)
        for j in range(self.M):
            if(self._fillingScheme[j]):
                self.A[j]=self.A_GLOBAL
                self._bunchLength[j]=self.BUNCH_LENGTH_GLOBAL
                self.phi[j]=self.PHI_GLOBAL
         
        self._setBunches()
 

#%%           
    
    def plotLongitudinalProfile(self, tmin=-1, tmax=-1, fontsize=12 ,legendsize=15):
        [t,s]=self.longitudinalProfile
        if tmin==-1:
            tmin=np.min(t)
        if tmax==-1:
            tmax=np.max(t)

        
        #        if (self._fillNumber!=0):
#            label='fill: '+str(self._fillNumber)
#        else:
        label=str(self.filledSlots) + ' bunches'
        
        
        plt.figure()
        plt.plot(t*1e6,s/1e6,label=label)
        plt.grid(True)
        plt.ylim(0,np.max(s)/1e6)
        plt.xlim(tmin*1e6, tmax*1e6)
            
        plt.tick_params(axis='both', which='major')
        plt.xlabel(r't $[\mu s]$')
        plt.ylabel(r'Longitudinal time distribution $[1/\mu s]$')
        plt.legend(fontsize=legendsize, loc='upper right')
        plt.show()
        
#%%       
        
    def setNbFromFillNumber(self):
        
        db=pytimber.LoggingDB()
        bunchIntensities='LHC.BCTFR.A6R4.B'+str(self._beamNumber)+':BUNCH_INTENSITY'
        if(self._fillNumber!=0):
            
            fill=db.getLHCFillData(self.fillNumber)
            
            for j in range(len(fill['beamModes'])):
                if(fill['beamModes'][j]['mode']==self.fillMode):
                    MODE=j
                
            t1=fill['beamModes'][MODE]['startTime']
            t2=fill['beamModes'][MODE]['endTime']
            
            
        A=db.get(bunchIntensities,t1,t2)
        timeStamps,A=A[bunchIntensities]
        i=len(A)-1
        while True:
            if np.sum(A[i])!=0:
                break
            i=i-1
        A=A[i]
        mask=A!=0
        self.Nb=np.mean(A[mask])
        self._NbIsComputed=True
        print ("Nb updated: " , self.Nb/1e11, "e11")
        print ("Nb calculated at: ", pytimber.dumpdate(t2))


        
#%%

    def plotPowerSpectrum(self, fmin=0, fmax=2):
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
#        plt.legend(fontsize=20)
        plt.show()

#%%         

    def _plotSpectrumFromTimber(startDate,beamNumber):
        db=pytimber.LoggingDB()
        
        if(beamNumber==1):
            freq='ALB.SR4.B1:SPEC_FREQ'
            powSpec='ALB.SR4.B1:SPEC_POW'
        elif(beamNumber==2):
            freq='ALB.SR4.B2:SPEC_FREQ'
            powSpec='ALB.SR4.B2:SPEC_POW'
        
        
    
        delayTime=5*60
        
    
        t1 = pytimber.parsedate(startDate) -  delayTime
        t2= 'next'
        #t2= 'next'#'2017-05-20 14:55:00.000'
    
    
        
        f=db.get(freq,t1,t2)  #Siccome i filledbuckets non cambiano prendo tutti i valori dall'inizio 
                                         #al FLATTOP e poi seleziono solo l'ultimo vettore
        timeStamps,f=f[freq]
        f=f[-1]/1e9
    
        S=db.get(powSpec,t1,t2)
        timeStamps,S=S[powSpec]
        print('PowerSpectrum Date :', pytimber.dumpdate(timeStamps[-1]))
        S=S[-1]
    #    pl.figure()
    #    pl.plot(f, S, color='b',label=pytimber.dumpdate(timeStamps[-1]))
    #    pl.xlim(0,2)
    #    pl.ylim(np.min(S),np.max(S))
    #    pl.grid('on')
    #    pl.tick_params(axis='both', which='major', labelsize=14)
    #    pl.xlabel('f [GHz]', fontsize=18)
    #    pl.ylabel('Power spectrum [a.u.]', fontsize=18)
    #    pl.legend(fontsize=14)
    #    pl.savefig('Figures/NoisePowerSpectrumSemiLog.png', dpi=200)
    #    pl.show()
        
        
#    
        return f,S,timeStamps[-1]
        
#%%  

    def getPloss(self,Z):
#        if(self._NbIsComputed==False):
#            self.setNbFromFillNumber()
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
            
            
#            pl.figure()
#            pl.plot(f[mask]/1e9,S,label='normalized spectrum')
#            pl.plot(f[mask]/1e9,Z.real, color='r', label='Impedance [ohm]')
#            pl.xlim(0,2)
#            pl.ylim(0,)
#            pl.tick_params(axis='both', which='major', labelsize=14)
#            pl.xlabel("f [GHz]",fontsize=18)
#            pl.ylabel("Impedance",fontsize=18)
#            pl.grid('on')
#            pl.legend()
#            pl.show()
    
        mask3= f>=0
        f=f[mask3]
        S=S[mask3]   
        Zreal=Zreal[mask3]
        Zf=f

#        print 'I0: ',(f0*e*self.filledSlots*self.Nb)**2
        P=f0*e*S*self.filledSlots*self.Nb
        
        P_density=2*(P**2)*Zreal
        P_loss=np.sum(P_density)  
#        print Ploss
        return P_loss, P_density;


    def get_2_beam_Ploss(self, Z, phase_shift):
#        if(self._NbIsComputed==False):
#            self.setNbFromFillNumber()
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
            
            
#            pl.figure()
#            pl.plot(f[mask]/1e9,S,label='normalized spectrum')
#            pl.plot(f[mask]/1e9,Z.real, color='r', label='Impedance [ohm]')
#            pl.xlim(0,2)
#            pl.ylim(0,)
#            pl.tick_params(axis='both', which='major', labelsize=14)
#            pl.xlabel("f [GHz]",fontsize=18)
#            pl.ylabel("Impedance",fontsize=18)
#            pl.grid('on')
#            pl.legend()
#            pl.show()
    
        mask3= f>=0
        f=f[mask3]
        S=S[mask3]   
        Zreal=Zreal[mask3]
        Zf=f

#        print 'I0: ',(f0*e*self.filledSlots*self.Nb)**2
        P=f0*e*S*self.filledSlots*self.Nb
        P_loss = []
        P_density_list = []
        for shift in phase_shift:
            P_density=4*(P**2)*Zreal*(1-np.cos(2*np.pi*f*shift))
#            P_density_list.append(P_density)
            P_loss.append(np.sum(P_density))  
#        print Ploss
        return P_loss;# f, P_density_list;

#%%
 
    def setBeamsFromSumWithShift(self, beam1, beam2, shift):
        
        [t1,s1]=beam1.longitudinalProfile
        [t2,s2]=beam2.longitudinalProfile
        deltaT=t1[1]-t1[0]
        step=int(shift/deltaT)
        s1=np.roll(s1,step)
        s2=np.roll(s2,-step)
        
        self.longitudinalProfile = [t1,(s1+s2)/2]
        self.filledSlots = beam1.filledSlots + beam2.filledSlots
#        self.plotLongitudinalProfile()
        
        plt.figure()
        plt.plot(t1*1e6,s1,label='beam1',color='b')
        plt.plot(t2*1e6,s2,label='beam2',color='r',alpha=0.7)
#        plt.plot(t2*1e6,(s2+s1)[mask],label='beam1+beam2',color='g',alpha=0.7)
        plt.grid('on')
        plt.ylim(0,)
        plt.xlim(np.min(t1)*1e6,np.max(t1)*1e6)
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.xlabel("t [us]",fontsize=18)
        plt.ylabel("Longitudinal Profile [a.u.]",fontsize=18)
        plt.legend(fontsize=16)
        plt.show()
        
        
#%%
    def saveLongitudinalDistribution(self, path, phase_shift=0, normalise=False, plot=False):
        '''phase_shift is intendet to be in ns '''
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
        if plot:
            label=str(self.filledSlots) + ' bunches'
            fontsize=12
            legendsize=15
            plt.figure()
            plt.plot(t*1e6,s,label=label,color='b')
            plt.grid(True)
            plt.ylim(0,np.max(s))
            plt.xlim(np.min(t)*1e6, np.max(t)*1e6)

            plt.tick_params(axis='both', which='major', labelsize=fontsize)
            plt.xlabel(r't $[\mu s]$',fontsize=fontsize)
            plt.ylabel(r'Longitudinal time distribution $[1/\mu s]$',fontsize=fontsize)
            plt.legend(fontsize=legendsize, loc='upper right')
            plt.show()

#%% 

    def __collide(beam1,beam2):
        
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
#            plt.plot(t2[mask]*1e6,(s2+s1)[mask],label='beam1+beam2',color='g',alpha=0.7)
            plt.grid('on')
            plt.ylim(0,)
            plt.xlim(np.min(t1)*1e6,np.max(t1)*1e6)
            plt.tick_params(axis='both', which='major', labelsize=14)
            plt.xlabel("t [us]",fontsize=18)
            plt.ylabel("Longitudinal Profile [a.u.]",fontsize=18)
            plt.legend(fontsize=16)
            plt.show()
        
        
        return ;


#%%
###############################################################################
###############################################################################

class Impedance(object):
    
    def __init__(self, f, CST_file=0):
        self.isResonatorImpedance=False
        self.isRWImpedance=False
        self.f=f
        self.df = pd.DataFrame()
        self.df['f']=f
        self.df['Zr']=np.zeros(len(f))
        if CST_file!=0:
            self.getImpedanceFromCST(CST_file)
        
#%%        
        
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

#%%        
        
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

#%%
        
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
    

#%%    
    
    def plot(self,fMin=0,fMax=2e9, fontsize=12, legendsize=15):
        
        
        plt.figure()
        plt.plot(self.df['f']/1e9,self.df['Zr'],label='Re[Z]')
#        pl.plot(f,self.imag,color='r',label='Im[Z]',alpha=0.7)
        plt.tick_params(axis='both', which='major')
        plt.xlim(fMin/1e9,fMax/1e9)
        plt.ylim(0, )
        plt.title('Impedance')
        plt.xlabel('f [GHz]')
        plt.ylabel(r'Z $[\Omega]$')
        plt.grid(True)
        plt.legend()
        plt.show()
    
#
#
#[s,t]=beam1.longitudinalProfile
#
#beam2=Beam()
#beam2.setBeamFromFillNumber(5979,'FLATTOP',1)
#beam2.plotLongitudinalProfile()
#beam1.plotLongitudinalProfile()
#
#beam=Beam()
#beam.setBeamFromFillNumber(5979,'FLATTOP',1)
#Z=Impedance()
#[S,f]=beam.getSpectrum()
#Z.getResonatorImpedance(f,1,1e3,500e6)
#

