import numpy as np
import pytimber
import sys
import math

from bihc.impedance import Impedance
from bihc.power import Power
from bihc.plot import Plot

class Beam(Impedance, Power, Plot):
    '''Defining the beam characteristics.

    Enables to create a beam with certain characteristics. The beam can be
    created from an existing fill stored in Timber or from a self defined
    filling scheme. It is possible to assign the bunch lenght, the bunch
    shape, the total intensity.
    '''

    def __init__(self, M=3564, A=1, fillNumber=0,
                 bunchLength=1.2e-9, phi=0, realMachineLength=True,
                 ppbk=250,d=0 , Nb=2.3e11, bunchShape='GAUSSIAN', qvalue=1.2, beamNumber=1, 
                 fillMode='FLATTOP', fillingScheme=[False]*3564, machine='LHC'):
        
        c = 299792458 # Speed of light in vacuum [m/s]
        
        self.M = M # Default max numebr of buckets
        self.A_GLOBAL = A
        self.BUNCH_LENGTH_GLOBAL = bunchLength/4 # Bunch lenght [s]
        
        self.PHI_GLOBAL = phi
        self.realMachineLength = realMachineLength
        self.ppbk = ppbk

        self._bunchShape = bunchShape # Bunche shape (analytical function)
        self._q = qvalue #q value for the q-gaussian distribution 1>q>3
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
   
    @property
    def bunchLength(self):
        return self._bunchLength
    
    @bunchLength.setter
    def bunchLength(self,newBunchLength):
        print('updating bunch length')
        self._bunchLength=np.zeros(self.M)
        self._bunchLength[self._fillingScheme]=newBunchLength/4
        self._setBunches()

    @property
    def bunchShape(self):
        return self._bunchShape
    
    @bunchShape.setter
    def bunchShape(self,newShape):
        shapeList=['GAUSSIAN', 'BINOMIAL' , 'COS2', 'q-GAUSSIAN']
        print('updating bunch shape...')
        if newShape in shapeList:
            self._bunchShape=newShape
            self._setBunches()
            self._isSpectrumReady=False
        else:
            raise ValueError("bunchShape should be: 'GAUSSIAN', 'BINOMIAL', 'q-GAUSSIAN', or 'COS2'." )

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

    @property
    def spectrum(self):
        if self._isSpectrumReady:
            [f.s] = self._spectrum
            self.powerSpectrum=[f, np.abs(s)**2]
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
            self.powerSpectrum=[f, np.abs(S)**2]
            self._isSpectrumReady=True
            
            return self._spectrum
        
    @spectrum.setter
    def spectrum(self, newSpectrum):
#        self._spectrum=newSpectrum
        raise Exception("Spectrum can not be assigned")
                       
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

                    #set to zero the part of each bunch that is outside the l range arount his mean
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
                    
                    def _Cq(q):
                        Gamma = math.gamma
                        if q<=0.99415629720:
                            return (2.0*np.sqrt(np.pi))/((3.0-q)*np.sqrt(1-q))*(Gamma(1.0/(1.0-q)))/Gamma((3.0-q)/2.0/(1.0-q))
                        elif (q<1.005827 and q>0.99415629720):
                            return np.sqrt(np.pi)
                        elif (q>=1.005827 and q<3.0):
                            return (np.sqrt(np.pi)*Gamma((3.0-q)/2.0/(q-1.0)))/(np.sqrt(q-1.0)*Gamma(1.0/(q-1.0)))
                        else:
                            raise Exception('q>3.0')
                     
                    def _eq(x,q):
                        eq = np.zeros(len(x))
                        for i,xx in enumerate(x):
                            if ((q!=1) and (1+(1-q)*xx)>0):
                                eq[i] = (1+(1-q)*xx)**(1/(1-q))
                            elif q==1:
                                eq[i] = np.exp(xx)
                            else:
                                eq[i] = 0
                        return eq

                    #q-Gaussian function
                    b = self._bunchLength[i]
                    q = self.q #[TODO] pass 1 q-val per bunch
                    mu = self.phi[i]
                    sTemp=np.sqrt(b)/_Cq(q)*_eq(-b*(tTemp-mu)**2,q)  

                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(sTemp>0)
                    sTemp=sTemp*mask2*mask
                    
            else:
                sTemp=np.zeros(len(tTemp))
            if(i==0):
                s=sTemp
            else:
                s=np.concatenate((s,sTemp),axis=0)
        
        
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
        for j in range(len(fill['beamModes'])):
            if(fill['beamModes'][j]['mode']==fillMode):
                MODE=j
            
        t1=fill['beamModes'][MODE]['startTime']
        t2=fill['beamModes'][MODE]['endTime']
        print ('Mode selected: ',self.fillMode, 'starts at: ', pytimber.dumpdate(t1), 'ends at',pytimber.dumpdate(t2))
    
        fb=db.get(filledBuckets,ts,t2)
        timeStamps,fb=fb[filledBuckets]
  
        std=db.get(bunchLengths,t1,t2)
        timeStamps,std=std[bunchLengths]
        
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
        
    def setCustomBeam(self):

        self.isATimberFill=False

        self._fillingScheme[0:self.M]=[True]*self.M
        self._bunchLength[self._fillingScheme]=self.BUNCH_LENGTH_GLOBAL  #std vector of a single turn in the machine
        self.phi=np.ones(self.M)*self.PHI_GLOBAL
            
        self._setBunches()

        
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
 
    def setBeamsFromSumWithShift(self, beam1, beam2, shift):
        
        [t1,s1]=beam1.longitudinalProfile
        [t2,s2]=beam2.longitudinalProfile
        deltaT=t1[1]-t1[0]
        step=int(shift/deltaT)
        s1=np.roll(s1,step)
        s2=np.roll(s2,-step)
        
        self.longitudinalProfile = [t1,(s1+s2)/2]
        self.filledSlots = beam1.filledSlots + beam2.filledSlots
        

        