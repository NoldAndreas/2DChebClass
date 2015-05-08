%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

nParticlesS=[20;20;20];

kBT=1;          % temperature
mS=[1;1;1];

gammaS=[1;1;1];
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

%V1DV1='free2D_Box';
%V1DV1 = 'quad_Box';

% appropriate physical parameters for potentials in V1DV1
% L1S = 8;
% L2S = 6;
% kBTS = kBT;
% 
% tauS = 1;
% y10S = 5;
% y20S = 5;
% B10S = 1;
% B20S = 1;
% V0AddS  = 1;
% potParamsNames = {'L1','L2', ...
%                   'tau','y10','y20','B10','B20','V0Add'};

%V1DV1='V1_Test_Box_Eq';
V1DV1='V1_Test_Box_Move';

% appropriate physical parameters for potentials in V1DV1

y0 = 10;

% box
L1S = y0;
L2S = y0;
kBTS = kBT;

% potential
V0S = 2;
y10S = 8;
y20S = 7;
tauS = 10;


% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'L1','L2', ...
                  'tau','y10','y20','V0'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='Gaussian';

 epsilonS=2*[ 1 1 1 ;
                 1 1 1 ;
                 1 1 1 ] ;

alpha1 = 0.5;
alpha2 = 1;
alpha3 = 1.5;

alpha12=(alpha1+alpha2)/2;
alpha13=(alpha1+alpha3)/2;
alpha23=(alpha2+alpha3)/2;

alphaS=[alpha1  alpha12 alpha13 ; 
        alpha12 alpha2  alpha23 ;
        alpha13 alpha23 alpha3  ];


potParams2Names={'epsilon','alpha'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax = 5;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true
%nSamples=100000;  

nSamples=500000;  

sampleFinal = false;

initialGuess='makeGridPosScale';

% number of runs of stochastic dynamics to do, and average over
%nRuns=50000;
nRuns=100;

% number of cores to use in parallel processing
poolsize=12;
%poolsize=4;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r','rv','r','rv'};

% whether to include hydrodynamic interactions
stocHI={false,false,true,true};
% HI interaction matrices
stocHIType={[],[],'RP','RPInv'};

% names for stochastic calculations -- used as legend text
stocName={'r0','rv0','r1','rv1'};

% whether to do Langevin and Brownian dynamics
%doStoc={true,true,true,true};
doStoc={false,false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true,true};

% number of time steps
tSteps={10^4,10^3,2*10^4,10^3};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};

stocColour = {{'r','b','g'},{'g'},{'b'},{'m'}};
stocStyle = {{'--','--','--'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

y0 = 3;


Phys_Area = struct('shape','Box','L1',L1S,'L2',L2S);

Plot_Area = struct('y1Min',0,'y1Max',L1S,'N1',100,...
                       'y2Min',0,'y2Max',L2S,'N2',100);

V2_Num   = struct('Fex','Meanfield','N',[20;20],'L',1);

paramNames = {'PhysArea','PlotArea','V2Num','doPlots'};

PhysArea = {};
PlotArea = {};
V2Num = {};
DDFTCode = {};
DDFTParamsNames = {};
DDFTType = {};


DDFTName = {};
doDDFT = {};
loadDDFT = {};

for N = 10:2:50
    Phys_Area.N = [N,N];
    PhysArea = cat(2,PhysArea,Phys_Area);
    PlotArea = cat(2,PlotArea,Plot_Area);
  %  FexNum = cat(2,FexNum,Fex_Num);
    V2Num = cat(2,V2Num,V2_Num);
    DDFTCode = cat(2,DDFTCode,'DDFTDynamics');
    DDFTParamsNames = cat(2,DDFTParamsNames,{paramNames});
    DDFTType = cat(2,DDFTType,'r');
    doDDFT = cat(2,doDDFT,true);
    loadDDFT = cat(2,loadDDFT,true);
    DDFTName = cat(2,DDFTName, num2str(N));
end

                   
        
doPlots = false;

% for stochastic comparison
%DDFTColour = {{'k','c','m'},{'k','c','m'},{'k','c','m'},{'k','c','m'}};
%DDFTStyle = {{':',':',':'},{'-.','-.','-.'},{'--','--','--'},{'-','-','-'}};

% for N computations
%DDFTColour = {{'r','r','r'},{'g','g','g'},{'b','b','b'},{'m','m','m'}};
%DDFTStyle = {{':',':',':'},{'-.','-.','-.'},{'-','-','-'},{'--','--','--'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'contour';

% x axis for position and velocity plots

%for stochastic comparison
%rMin=[0;0];
%rMax=[L1S;L2S];

% for N computation
rMin = [4;4];
rMax = [10;10];

pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=6;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[0;0];
RMMax=[max(L1S,L2S);max(L1S,L2S)];
PMMin=[-1;-1];
PMMax=[1;1];

% number of bins for histograming of stochastic data
%nBins=[100;100];
nBins=[50;50];

separateSpecies = true;

% determine which movies/plots to make
% distribution movies/plots
doMovieGif     = false;          % .gif movie
doInitialFinal = false;
doMeans        = false;
doEquilibria   = true;

%sendEmail = true;
