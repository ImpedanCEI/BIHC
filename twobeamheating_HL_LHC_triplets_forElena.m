%% Two beam induced heating

clear all
close all

nsig=2.5;   %troncamento gaussiana
G=sqrt(5);
F=1.2413;
step=2;
cd('D:\data\repository_extraction\LHC\x_Carlo')
load('orbit_LHC_collision_HL_LHC.mat')
Zp=[];
% for m=31:1:31
% if m==1;
% load('Fill2736_t_0.800h_450.0GeV_B1.mat')
% elseif m==2;
% load('Fill2736_t_0.800h_450.0GeV_B2.mat')
% elseif m==3;
% load('Fill2736_t_1.330h_4000.0GeV_B1.mat')
% elseif m==4;
% load('Fill2736_t_1.330h_4000.0GeV_B2.mat')
% elseif m==5;
% load('Fill2736_t_2.550h_4000.0GeV_B1.mat')
% elseif m==6;
% load('Fill2736_t_2.550h_4000.0GeV_B2.mat')
% elseif m==7;
% load('Fill2736_t_6.662h_4000.0GeV_B1.mat')
% elseif m==8;
% load('Fill2736_t_6.662h_4000.0GeV_B2.mat')
% elseif m==9;
% load('Fill2736_t_10.775h_4000.0GeV_B1.mat')
% elseif m==10;
% load('Fill2736_t_10.775h_4000.0GeV_B2.mat')
% elseif m==11;
% load('Fill2736_t_14.887h_4000.0GeV_B1.mat')
% elseif m==12;
% load('Fill2736_t_14.887h_4000.0GeV_B2.mat')
% elseif m==13;
% load('Fill2736_t_19.000h_4000.0GeV_B1.mat')
% elseif m==14;
% load('Fill2736_t_19.000h_4000.0GeV_B2.mat')
% 
% elseif m==15;
% load('Fill3000_t_0.660h_450.0GeV_B1.mat')
% elseif m==16;
% load('Fill3000_t_0.660h_450.0GeV_B2.mat')
% elseif m==17;
% load('Fill3000_t_1.230h_4000.0GeV_B1.mat')
% elseif m==18;
% load('Fill3000_t_1.230h_4000.0GeV_B2.mat')
% elseif m==19;
% load('Fill3000_t_1.880h_4000.0GeV_B1.mat')
% elseif m==20;
% load('Fill3000_t_1.880h_4000.0GeV_B2.mat')
% elseif m==21;
% load('Fill3000_t_2.553h_4000.0GeV_B1.mat')
% elseif m==22;
% load('Fill3000_t_2.553h_4000.0GeV_B2.mat')
% elseif m==23;
% load('Fill3000_t_3.227h_4000.0GeV_B1.mat')
% elseif m==24;
% load('Fill3000_t_3.227h_4000.0GeV_B2.mat')
% elseif m==25;
% load('Fill3000_t_3.900h_4000.0GeV_B1.mat')
% elseif m==26;
% load('Fill3000_t_3.900h_4000.0GeV_B2.mat')
% 
% elseif m==27;
% load('Fill3286_t_2.470h_450.0GeV_B1.mat')
% elseif m==28;
% load('Fill3286_t_2.470h_450.0GeV_B2.mat')
% elseif m==29;
% load('Fill3286_t_3.030h_4000.0GeV_B1.mat')
% elseif m==30;
% load('Fill3286_t_3.030h_4000.0GeV_B2.mat')
% elseif m==31;
% load('Fill3286_t_4.000h_4000.0GeV_B1.mat')
% elseif m==32;
% load('Fill3286_t_4.000h_4000.0GeV_B2.mat')
% elseif m==33;
% load('Fill3286_t_6.375h_4000.0GeV_B1.mat')
% elseif m==34;
% load('Fill3286_t_6.375h_4000.0GeV_B2.mat')
% elseif m==35;
% load('Fill3286_t_8.750h_4000.0GeV_B1.mat')
% elseif m==36;
% load('Fill3286_t_8.750h_4000.0GeV_B2.mat')
% elseif m==37;
% load('Fill3286_t_11.125h_4000.0GeV_B1.mat')
% elseif m==38;
% load('Fill3286_t_11.125h_4000.0GeV_B2.mat')
% elseif m==39;
% load('Fill3286_t_13.500h_4000.0GeV_B1.mat')
% elseif m==40;
% load('Fill3286_t_13.500h_4000.0GeV_B2.mat')
% end
% i=m
%% Construction of the beam distribution

ninj=10;
nslots=3564;
ntrain=4;
nbunches=72;
nbunches2=36;
nbunches3=18;
nbunch1=72*ntrain*ninj;
nbunch2=36*ntrain*ninj;
nbunch3=18*ntrain*ninj;
t0 = 25e-9;
batchS=7; %batch spacing in 25 ns slots
injspacing=37; %injection spacing in 25 ns slots

% spacing = input('Insert the bunch spacing in s[25e-9]: ');
% if isempty(spacing)
%     spacing = 1;
% end
%25 ns beam
PWLall=[];
PWLserall=[];
BS=200; %batch spacing in ns
% create 25 ns beam    
Np=1.2e11;
bt=ones(1,nbunches);
st=zeros(1,batchS);
stt=zeros(1,injspacing);
sc=zeros(1,nslots-(ntrain*nbunches*ninj+((ntrain-1)*(batchS)*ninj)+((1)*injspacing*(ninj))));
an1=[bt st bt st bt st bt stt];
an=[repmat(an1,1,ninj) sc];
nbunch=nbunch1;

%% reading the simulated impedance
%MKE=dlmread('impedance_a_13cm.txt','',2,0);
MKE=dlmread('Z_long_real_beamscreen.txt','',2,0);
MKEser=1.2;
%% Calculating the Fourier transform on the LHC beam

% nplus=nslots*ngiri;
n=[1:1:nslots];
frev=11.245e3;
wrev=2*pi*frev;

%fit of the impedance data
pstop=150000;
pstep=1;
xp=[pstep*frev:frev*pstep:pstop*frev];
Z=real(MKE);
Zser=(MKEser);

Zint=interp1(MKE(:,1),Z,xp);
Zintser=Zint.*Zser;

e=1.602e-19;
c=3e8;
Q=Np*e*nbunch;
%Q=e*sum(an);
q=Np*e;

sigma=0.075;
%sigma=mean(nonzeros(sigmaz));

%striplet=[13353:step:13363];
striplet=[13363:step:13413];
 
%Vectors inizialization
PWLeq=[];
PWLsereq=[];
Ivect=[];
lambdavect=[];
FFTbeamp=[];
lambdap=[];
pp=[];
PL3=[];
PL2=[];
PL1=[];

%Yb1int=interp1(S_b1(10410:4:10490),Y_b1(10410:4:10490),striplet);
%Yb2int=interp1(S_b2(10410:4:10490),Y_b2(10410:4:10490),striplet);

%Orbit in the common chamber
Yb1int=interp1(S_b1(10460:4:10620),Y_b1(10460:4:10620),striplet);
Yb2int=interp1(S_b2(10460:4:10620),Y_b2(10460:4:10620),striplet);



strip=13363:step:13413;
    PWL=[];
    PWL2=[];
PWLser=[];
    istrip=int32((strip-13363))/step+1;
    ts=2*strip/c;
    Yb1intfor=Yb1int(istrip)+Yb2int(istrip);
    Yb2intfor=Yb1int(istrip)-Yb2int(istrip);
for p=1:pstep:pstop
 
    percint=p/pstop*100%calculation advancement in %
    
  
   lambda=exp(-(p.*p.*wrev.*wrev.*sigma.*sigma)/(2.*c.*c));
FFTbeam=(1/(sum(an)))*lambda.*sum(an.*exp(1i.*n*p*wrev*t0));
PWlossr1=2*(Q^2)*frev*frev*sum((abs(FFTbeam).^2).*Zint(p/pstep));
PWlossr2=-2*(Q^2)*frev*frev*sum((abs(FFTbeam).^2).*Zint(p/pstep)*diag(cos(p*wrev*ts)));
PWlossr3=2*(Q^2)*frev*frev*sum((abs(FFTbeam).^2)*25.*Zint(p/pstep)*diag(Yb1intfor)*(diag(ones(1,length(ts)))-diag(cos(p*wrev*ts))));
PWloss=2*(PWlossr1+PWlossr2+PWlossr3); %total two beam power loss
PWloss2=2*(PWlossr1+PWlossr2);
PWlossser=2*PWlossr1; %neglecting coupling between beams
% lambdavect=[lambdavect lambda];
%lambdap=[lambdap lambda];
%FFTbeamp=[FFTbeamp FFTbeam]; 
%pp=[pp p*frev];
PWL=[PWL PWloss];
PWL2=[PWL2 PWloss2];
PWLser=[PWLser PWlossser];
end
for i=1:1:length(strip)
ind=[i:length(strip):length(PWL)]
    PL=sum(PWL(ind));
    PL2v=sum(PWL2(ind));

PL3=[PL3 PL];
PL2=[PL2 PL2v];
%PL1=[PL1 PLser];
end
PLser=sum(PWLser)*length(strip);
PL1=PLser;
%clear PWL PWLser lambda FFTbeam


%[X0,Y0] = intersections(sigma,S21plot,sigma,S21mobtt,1);
%Zp=[Zp PLser];
%%
%ltriplet=13363-13353;
ltriplet=13413-13363;
IntPWLser=((step*trapz(PL1)))*1.2/ltriplet*1782/1440;
IntPWL=((step*trapz(PL3)))*1.2/ltriplet*1782/1440;
IntPWL2=((step*trapz(PL2)))*1.2/ltriplet*1782/1440;

%% post processing for the case under study

weld49=[1.27 1.26 1.235 1.2 1.18 1.16];
weld59=[1.23 1.22 1.200 1.175 1.155 1.14];

rhov=[0.047 0.053 0.070 0.103 0.147 0.200]*1e-8;
IntPWLn=IntPWL*weld49./1.2.*sqrt(rhov)./sqrt(rho)
IntPWL2n=IntPWL2*weld49./1.2.*sqrt(rhov)./sqrt(rho)
IntPWLsern=IntPWLser*weld49./1.2.*sqrt(rhov)./sqrt(rho)

