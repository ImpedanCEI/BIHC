'''
Beam module to manage beam object creation
from Timber database or specified by custom 
filling scheme defined by the user.

The created beam consists in bunches allocated in 
the buckets specified by the boolean filling scheme
defined by user. For this, the longitudinal beam profile 
in time and spectrum and power spectrum are calculated
provided a beam shape and bunch length.

* date: 12/12/2022
* author: Francesco Giordano, Elena de la Fuente, Leonardo Sito
'''

import numpy as np
import sys
import math

from bihc.impedance import Impedance
from bihc.power import Power
from bihc.plot import Plot

class Beam(Impedance, Power, Plot):
    '''Defining the beam characteristics.

    Enables to create a beam object with certain characteristics. The beam can be
    created from an existing fill stored in Timber or from a self defined
    filling scheme. It is possible to assign the bunch lenght, the bunch
    profile shape, the total intensity and the number of beams (1 or 2).

    Parameters
    ----------
    M : int, default 3564
        Maximum number of bucket
    A : int, default 1
        Normalized amplitude for bunch profiles
    bunchLength : float, default 1.2e-9
        Beam total longitudinal bunch length in seconds (4*sigma) [s]
    bunchShape : str, default 'GAUSSIAN'
        Beam profile shape : 'GAUSSIAN', 'BINOMIAL' , 'COS2' or 'q-GAUSSIAN'         
    qvalue : float, default 1.2
        q-Gaussian q-value for the 'q-GAUSSIAN' beam profile opt
    phi : float, default 0
        Offset of the bunch profile distribution in time [s]
    t0 : float, default 0
        The time length (space) of one bucket
    Np : float, default 2.3e11
        Beam intensity in number of protons per bunch
    beamNumber : int, default 1
        Number of beams for the power loss computation (1 or 2)
    fillNumber : int, default 0
        Fill number relative to a particular beam fill of the machine
    fillMode : str, default 'FLATTOP'
        Timber label to extract data at a certain energy 'INJ', 'FLATTOP', 'STABLE'
    fillingScheme : list of bool, default [False]*3564
        Bool values to define the bunch filling scheme with length the number of buckets
    machine : str, default 'LHC'
        Name of the machine to operate with : 'PS', 'SPS', 'LHC'
    spectrum : str, default 'numeric'
        Whether to calculate the spectrum with a numerical FFT 'numeric', from the analytical formula 'analytic', or from input 'user'
    realMachineLength : bool, default True
        Flag to adapt bucket size to real machine length
    ppbk : int, default 250
        Number of time samples per bucket
    frev : float, default None
        Revolution frequency in [Hz] to sample the analytic beam spectrum computation
    fmax : float, default 2e9
        Maximum frequency in [Hz] up to which to compute the beam spectrum 
    verbose : bool, default True
        Flag to control console output

    Attributes
    ----------
    longitudinalProfile : numpy.ndarray list
        Beam longitudinal distribution in time, normalized. Returns the list of numpy arrays [time, profile]
    spectrum : numpy.ndarray list 
        Beam spectrum in frequency. Returns the list of numpy arrays [frequency, spectrum]
    powerSpectrum : numpy.ndarray list
        Beam power spectrum in frequency. Returns the list of numpy arrays [frequency, powerspectrum]
    totalBeamCharge : float
        Beam charge computed from intensity and number of filled slots, in Coulombs [C]
    '''

    def __init__(self, M=3564, A=1, fillNumber=0,
                 bunchLength=1.2e-9, phi=0, realMachineLength=True,
                 ppbk=250, t0=None, Np=2.3e11, bunchShape='GAUSSIAN', LPCfile=None, qvalue=1.2, beamNumber=1, 
                 fillMode='FLATTOP', fillingScheme=[False]*3564, machine='LHC', spectrum='numeric', frev=None, 
                 fmax=2e9, exp=2.5, verbose=False):
        
        c = 299792458 # Speed of light in vacuum [m/s]
        
        self.M = M # Default max numebr of buckets
        self.A_GLOBAL = A
        self.BUNCH_LENGTH_GLOBAL = bunchLength/4 # Bunch lenght (sigma) [s]
        
        self.PHI_GLOBAL = phi
        self.realMachineLength = realMachineLength
        self.ppbk = ppbk

        self._bunchShape = bunchShape # Bunche shape (analytical function)
        self.q = qvalue #q value for the q-gaussian distribution 1>q>3
        self.exp = exp
        self.J = 1
        self.Np = Np # Number of particles per bunch
        self.fillMode = fillMode
        self._fillNumber = fillNumber # Fill number relative to a particular fill of the machine
        
        self._isSpectrumReady = False
        self.isATimberFill    = False
        self.verbose = verbose 
        self._machine = machine # Select the machine you are working with
        self._spectrumtype = spectrum
        self.frev = frev
        self.fmax = fmax

        # Some parameters are set in a different way depending on the machine you are working with
        if self._machine == 'LHC':
            self.BUCKET_MAX = 3564
            RING_CIRCUMFERENCE = 26658.883   #[m]
            GAMMA_R = 7461                   # flat top
                        
        elif self._machine =='SPS':
            self.BUCKET_MAX = 924 
            RING_CIRCUMFERENCE = 6911          #[m]
            if self.fillMode == 'FLATTOP':
                GAMMA_R = 251    # flat top 450 GeV
            else:
                GAMMA_R = 27.7   # flat bottom value 26 GeV

        elif self._machine == 'PS':
            self.BUCKET_MAX = 21
            RING_CIRCUMFERENCE = 628 #2*pi*100
            GAMMA_R = 27.7366 #28.7185  # p=26GeV

        elif self._machine== 'PSB': #TODO
            pass
            
        if self.M > self.BUCKET_MAX:
            self.M = self.BUCKET_MAX
            if self.verbose:
                print('Number of bucket that could be filled set to: ', self.M)
        
        BETA_R = np.sqrt(1 - (1/GAMMA_R**2))
        self.T_1_TURN = RING_CIRCUMFERENCE/(c*BETA_R)
        if self.frev is None:
            self.frev = 1/self.T_1_TURN 

        # t0 is the time (space) of one bucket 
        if t0 is None:
            self.t0 = self.T_1_TURN/self.BUCKET_MAX #forced to be integer
        else: 
            self.t0 = t0
        
        self.l = self.t0/2
        self._bunchLength = np.zeros(self.M)
        self.phi = np.zeros(self.M)
        self._spectrum = [np.zeros(self.M), np.zeros(self.M)]    #[f,S]
        self._fillingScheme=fillingScheme[0:self.M]
             
        self._beamNumber = beamNumber # Beam number, either 1 or 2
        self._beamFile = LPCfile

        # Computes the beam longitudinal profile 
        if fillNumber > 0:
            #if user specifies a fill number, data is extracted from timber
            self.setBeamFromFillNumber(fillNumber, fillMode, beamNumber)
        elif sum(self._fillingScheme) > 0:
            #if the array is not all False values
            self.setCustomBeamWithFillingScheme()
        elif LPCfile is not None:
            self.LPCfile = LPCfile.split('.')[0]
            self.setBeamFromLPC()
        else:
            #if the array is all False values
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
        shapeList=['GAUSSIAN', 'BINOMIAL' ,'PARABOLIC', 'COS2', 'q-GAUSSIAN']
        print('updating bunch shape...')
        if newShape in shapeList:
            self._bunchShape=newShape
            self._setBunches()
            self._isSpectrumReady=False
        else:
            raise ValueError("bunchShape should be: 'GAUSSIAN', 'BINOMIAL', 'PARABOLIC', 'q-GAUSSIAN', or 'COS2'." )

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
        '''spectrum (property)

        Computes spectrum and power spectrum from the 
        longitudinal beam profile using Numpy fft or 
        the analytical fomula (C. Zannini)

        Returns
        -------
        spectrum : numpy.ndarray list 
            Beam spectrum in frequency. Returns the list of numpy arrays [frequency, spectrum]
        '''
        self.powerSpectrum = []
        if self._isSpectrumReady:
            [f,s] = self._spectrum
            self.powerSpectrum=[f, np.abs(s)**2]
            return self._spectrum

        elif self._spectrumtype == 'user':
            return self._spectrum #spectrum must be provided through the setter

        else:
    
            if self._spectrumtype == 'numeric':
                [t,s]=self.longitudinalProfile
                s=s/(self.filledSlots*self.J) #why is it normalized?
                deltaT=t[10]-t[9] 
                fc=1/deltaT                       #in frequency we have a periodic signal of period fc where fc is 1/deltaT where the sampling step
                                                  #We are interested only in the range between -fc/2 and fc/2 ( In particular (0,fc/2) beacouse x is real)
                S=np.fft.fft(s,len(s))
                S=np.fft.fftshift(S)
                S=S*deltaT
                deltaF=fc/len(s)
                f = np.linspace(-fc/2, fc/2 - deltaF, len(S))      #vector of K point from min_value to max_value (Domain in frequency)

            else: #analytic formula (C.Zannini)
                print('\033[93m'+'Warning -> Analytic FFT might take a long time to be computed'+'\033[0m')
                from  bihc.plot import progressbar

                an = self._fillingScheme
                t0 = self.t0             #25 ns
                n = np.arange(1,self.M+1)
                c = 299792458
                wrev = 2*np.pi*self.frev
                sigma = self.BUNCH_LENGTH_GLOBAL*c #sigma in m
                sigmacos= 0.854*sigma #match FWHM to the gaussian sigma
                sigmapar= 0.744653*sigma
                F = 1.2413 
                A = (1/np.sum(an))

                #S = np.zeros(self.M*self.ppbk) #same length as numeric
                S = np.zeros(int(self.fmax/self.frev))
                lambdas = np.zeros_like(S)

                if(self._bunchShape=='BINOMIAL'): #todo
                    raise Exception("BINOMIAL is not supported for analytic spectrum calculation")
                
                elif(self._bunchShape=='GAUSSIAN'):
                    for p in progressbar(range(len(S)), "Computing analytic FFT: ", 20):
                        lambdas[p]=np.exp(-(p*p*wrev*wrev*sigma*sigma)/(2*c*c))
                        S[p]=np.abs(A*lambdas[p]*np.sum(an*np.exp(1j*p*wrev*n*t0)))
            
                elif(self._bunchShape=='COS2'):
                    for p in progressbar(range(1,len(S)), "Computing analytic FFT: ", 20): #TODO fix
                        Fc=(F**2)*(sigmacos**2)*((p*wrev)**2)/(c**2)
                        lambdas[p]=-1.14*np.sqrt(2*np.pi)/np.pi/(sigmacos*p*wrev/c*(-2+Fc))*(np.sqrt(2/np.pi))*np.sin((np.pi*sigmacos*p*wrev*F)/(np.sqrt(2)*c))
                        S[p]=np.abs(A*lambdas[p]*np.sum(an*np.exp(1j*p*wrev*n*t0)))

                elif(self._bunchShape=='PARABOLIC'):
                    for p in progressbar(range(1,len(S)), "Computing analytic FFT: ", 20): 
                        CosSin=(np.sqrt(5)*sigmapar*p*wrev/c*np.cos(np.sqrt(5)*sigmapar*p*wrev/c)-np.sin(np.sqrt(5)*sigmapar*p*wrev/c))
                        lambdas[p]=-3*c**3/((np.sqrt(5)**3)*(sigmapar**3)*(p*wrev)**3)*CosSin
                        S[p]=np.abs(A*lambdas[p]*np.sum(an*np.exp(1j*p*wrev*n*t0)))
                
                f = np.linspace(1,len(S),len(S))*self.frev #[TODO] should it start in 0?
                self.lambdas=[f, lambdas]

            if self.verbose:
                print ('DC component: ', np.max(np.abs(S)))  

            self._spectrum=[f, np.abs(S)]
            self.powerSpectrum=[f, np.abs(S)**2]
            self._isSpectrumReady = True
            
            return self._spectrum
        
    #@spectrum.setter
    def setSpectrum(self, newSpectrum): 
        self._spectrum = newSpectrum
        self._isSpectrumReady = True
        #raise Exception("Spectrum can not be assigned")
                       
    def _setBunches(self):
        '''_setBunches method

        Computes the longitudinal profile of the bunches 
        with the shape specified in the class instance

         'GAUSSIAN', 'BINOMIAL' , 'COS2' or 'q-GAUSSIAN' 
        '''
        if((self.realMachineLength) and (self.t0*self.M>self.T_1_TURN)):
            self.t0=self.T_1_TURN/self.M
            if self.verbose:
                print("t0 has been resized to ", self.t0 ," because it was to big to fill ",self.M," buckets in the real machine length")
        
        deltaD=self.t0/self.ppbk
        tTemp = np.arange(-self.t0/2, self.t0/2 ,deltaD)
        
        print ("Elaborating Data...")
        
        s = np.zeros(len(tTemp)) #[0] * len(t)
        sTemp=np.zeros(len(s))
        H=0
        self.filledSlots=0
        
        for i in range(self.M):
            if self._fillingScheme[i]==True:
                
                self.filledSlots+=1 #TODO: how to match the bunch length
                if(self._bunchShape=='BINOMIAL'):
                    H=(np.sqrt(2/np.log(2)))*(self._bunchLength[i]*4) #Binomial function (Francesco)
                    #H=(np.sqrt(2/np.log(2)))*(self._bunchLength[i]*2*np.sqrt(2*np.log(2)))
                    sTemp=(1 - 4*((tTemp - self.phi[i])/(H))**2) #Binomial function (Francesco)

                    #set to zero the part of each bunch that is outside the l range arount his mean
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(sTemp>0)
                    sTemp=sTemp*mask2*mask
                    sTemp=2*(sTemp**self.exp)/H
                    sTemp=np.sum(sTemp)*deltaD
                    profile_1_bunch = [tTemp, sTemp]

                elif(self._bunchShape=='GAUSSIAN'):
                    sTemp=1/(self._bunchLength[i] * np.sqrt(2 * np.pi)) *(np.e)**(-((tTemp - self.phi[i])**2)/(2*self._bunchLength[i]**2))  #Gaussian function
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=np.ones(len(sTemp))
                    sTemp=sTemp*mask2*mask
                    profile_1_bunch = [tTemp, sTemp]
                    
                elif(self._bunchShape=='COS2'):
                    tc=self._bunchLength[i]*2.77 #4*0.854
                    sTemp=1/tc*(np.cos(np.pi*(tTemp - self.phi[i])/(2*tc)))**2  #cos^2 function
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(abs(tTemp)<tc)
                    sTemp=sTemp*mask2*mask
                    profile_1_bunch = [tTemp, sTemp]

                elif(self._bunchShape=='PARABOLIC'):
                    sTemp=(1-(1/(4*0.744653*self._bunchLength[i]**2))*(tTemp- self.phi[i])**2)  #Parabolic function 
                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(sTemp>0)       
                    sTemp=sTemp*mask*mask2
                    sTemp=sTemp/(sum(sTemp)*deltaD)
                    profile_1_bunch = [tTemp, sTemp]

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
                    q = self.q #[TODO] pass 1 q-val per bunch
                    b = 1/((self._bunchLength[i]**2)*(5-3*q))
                    mu = self.phi[i]
                    sTemp=np.sqrt(b)/_Cq(q)*_eq(-b*(tTemp-mu)**2,q)  

                    mask=(np.abs(tTemp - self.phi[i]))<self.l
                    mask2=(sTemp>0)
                    sTemp=sTemp*mask2*mask
                    profile_1_bunch = [tTemp, sTemp]
                    
            else:
                sTemp=np.zeros(len(tTemp))
            if(i==0):
                s=sTemp
            else:
                s=np.concatenate((s,sTemp),axis=0)
        
        
        t_max=self.t0*self.M-self.t0/2
        
        if(self.realMachineLength):
            t=np.arange(-self.t0/2,self.T_1_TURN -self.t0/2,deltaD)     
        else:
            t=np.arange(-self.t0/2,t_max,deltaD)
         
        s_end=np.zeros(len(t) - len(s))
        s=np.concatenate((s,s_end),axis=0)
        sTemp=s
        for k in range(1,self.J):
            s=np.concatenate((s,sTemp),axis=0)
        
        tn=t
        t=np.linspace(np.min(tn),np.max(tn),len(s))
        
        self.totalBeamCharge=self.Np*1.602176e-19*self.filledSlots
        #s=s/(self.filledSlots*self.J) #??? why is it normalized
        dt=t[1]-t[0]
    
        s_integr=np.sum(s)*dt
    
        mask=self._bunchLength!=0
        if self.verbose:
            print ('Average bunch length = ', np.mean(self._bunchLength[mask])*4, 's')
            print ('Integral of s = ' , s_integr)
            print ("Buckets filled : ", self.filledSlots)
            print ('Np: ',self.Np)
            print ("Total beam charge: ", self.totalBeamCharge, "C")
        
        self.profile_1_bunch = profile_1_bunch
        self.longitudinalProfile=[t,s]


    def setBeamFromFillNumber(self, fillNumber, fillMode='FLATTOP', beamNumber=1):
        '''Set beam from fill number

        Retrieves beam fill information from Timber provided a fill Number
        It requires to install `pytimber` python package
        For more information refer to bihc installation guide

        **Parameters will override class instantiation**

        Parameters
        ----------
        beamNumber : int, default 1
            Number of beams for the power loss computation (1 or 2)
        fillNumber : int, default 0
            Fill number relative to a particular beam fill of the machine
        fillMode : str, default 'FLATTOP'
            Timber label to extract data at a certain energy 'INJ', 'FLATTOP', 'STABLE'

        Raises
        ------
        Exception
            pytimber could not be imported
        '''
        try:
            import pytimber
        except:
            print('This method uses pytimber. Please follow the installation guide to set it in your python environment')

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
        if self.verbose:
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
        if self.verbose:
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
        self.setNpFromFillNumber()
        self._setBunches()
        
    def setCustomBeam(self):
        '''Set custom beam without a filling scheme

        Sets beam with all bunches set to True for all bucket 
        slots with the bunch length, bunch shape and offset 
        defined in class instance
        '''
        self.isATimberFill=False

        self._fillingScheme[0:self.M]=[True]*self.M
        self._bunchLength[self._fillingScheme]=self.BUNCH_LENGTH_GLOBAL  #std vector of a single turn in the machine
        self.phi=np.ones(self.M)*self.PHI_GLOBAL
            
        self._setBunches()

    def setBeamFromLPC(self):
        '''Set beam from LPC tool csv output

        Sets beam reading the rows of the .csv file 
        specified by the user. This .csv file is the
        ouput of the graphical tool LPC

        For more information check the tool documentation:
        https://lpc.web.cern.ch/schemeEditor.html
        '''
        import csv

        self.isATimberFill = False
        self._fillingScheme = np.zeros(self.M, dtype=bool)

        with open(self._beamFile) as f:
            data = csv.reader(f)
            for row in data:
                if row[0]: #avoid empty rows
                    self._fillingScheme[int((int(row[0])-1)/10)] = True

        self._bunchLength[self._fillingScheme]=self.BUNCH_LENGTH_GLOBAL  #std vector of a single turn in the machine
        self.phi=np.ones(self.M)*self.PHI_GLOBAL

        self._setBunches()

    def setCustomBeamWithFillingScheme(self):
        '''Set custom beam with a filling scheme

        Sets beam with bunches where the filling scheme list 
        is set to true, and empty when False; for every bucket
        slots, and with the bunch length, bunch shape and offset 
        defined in class instance
        '''
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
        
    def setNpFromFillNumber(self):
        '''Set Intensity from fill number

        Retrieves intensity  information from Timber provided a fill Number
        It requires to install `pytimber` python package
        For more information refer to bihc installation guide

        **Parameters will override class instantiation**

        Raises
        ------
        Exception
            pytimber could not be imported
        '''
        try:
            import pytimber
        except:
            print('This method uses pytimber. Please follow the installation guide to set it in your python environment')

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
        self.Np=np.mean(A[mask])
        self._NpIsComputed=True
        if self.verbose:
            print ("Np updated: " , self.Np/1e11, "e11")
            print ("Np calculated at: ", pytimber.dumpdate(t2))
 
    def setBeamsFromSumWithShift(self, beam1, beam2, shift):
        ''' Set beam object from the sum of two beam objects
        
        Parameters
        ----------
        beam1 : Beam object
            First beam to add
        beam2 : Beam object
            Second beam to add
        shift : float
            Time shift between beam 1 and beam 2 in seconds [s]
        '''
        [t1,s1]=beam1.longitudinalProfile
        [t2,s2]=beam2.longitudinalProfile
        deltaT=t1[1]-t1[0]
        step=int(shift/deltaT)
        s1=np.roll(s1,step)
        s2=np.roll(s2,-step)
        
        self.longitudinalProfile = [t1,(s1+s2)/2]
        self.filledSlots = beam1.filledSlots + beam2.filledSlots
        

        
