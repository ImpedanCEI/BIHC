# Define HL filling scheme
ninj=10
nslots=3564
ntrain=4
nbunches=72
nbunches2=36
nbunches3=18
nbunch1=72*ntrain*ninj
nbunch2=36*ntrain*ninj
nbunch3=18*ntrain*ninj
t0 = 25e-9;
batchS=7; # batch spacing in 25 ns slots
injspacing=37; # injection spacing in 25 ns slots

PWLall=[]
PWLserall=[]
BS=200 # batch spacing in ns

Np=1.2e11
bt=[True]*nbunches
st=[False]*batchS
stt=[False]*injspacing
sc=[False]*(nslots-(ntrain*nbunches*ninj+((ntrain-1)*(batchS)*ninj)+((1)*injspacing*(ninj))))
an1=bt+ st +bt+ st+ bt+ st+ bt+ stt
# an=[repmat(an1,1,ninj) sc] # This is the final true false sequence that is the beam distribution
an=an1 * ninj + sc

b_HL_2760b = Beam(bunchShape='GAUSSIAN', beamNumber=1, fillingScheme=an, bunchLength=1.2e-9, Nb=2.1e11, d=25e-9)
b_HL_2760b.plotLongitudinalProfile()
[f,S] = b_HL_2760b.spectrum    
b_HL_2760b.plotPowerSpectrum()
