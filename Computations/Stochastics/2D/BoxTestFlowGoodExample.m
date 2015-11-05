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

V1DV1='V1_Test_Box_Flow';

% appropriate physical parameters for potentials in V1DV1

y0 = 10;

% box
L1S = 30;
L2S = 5;
kBTS = kBT;

% potential
V0S = 0.1;
y10S = L1S/4;
y20S = L2S/2;
tauS = 1;

alpha1w = 1;
alpha2w = 1;
alpha3w = 1;

epsilon1w = 0.1;
epsilon2w = 2;
epsilon3w = 10;

alphaWallS = [alpha1w;alpha2w;alpha3w];
epsilon0S  = [epsilon1w;epsilon2w;epsilon3w];

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'L1','L2', ...
                  'tau','y10','y20','V0','alphaWall','epsilon0'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='Gaussian';

epsilon1 = 0.1;
epsilon2 = 0.1;
epsilon3 = 0.1;

epsilon12 = (epsilon1+epsilon2)/2;
epsilon13 = (epsilon1+epsilon3)/2;
epsilon23 = (epsilon2+epsilon3)/2;

epsilonS= [ epsilon1 epsilon12 epsilon13 ;
            epsilon12 epsilon2 epsilon23 ;
            epsilon13 epsilon23 epsilon3] ;

% alpha1 = 0.5;
% alpha2 = 1;
% alpha3 = 1.5;

alpha1 = 1;
alpha2 = 1;
alpha3 = 1;


alpha12=(alpha1+alpha2)/2;
alpha13=(alpha1+alpha3)/2;
alpha23=(alpha2+alpha3)/2;

alphaS=[alpha1  alpha12 alpha13 ; 
        alpha12 alpha2  alpha23 ;
        alpha13 alpha23 alpha3  ];


potParams2Names={'epsilon','alpha'};

%--------------------------------------------------------------------------
% U parameters
%--------------------------------------------------------------------------

UDU = 'Flow_Test_Box';

BoxHeightS = L2S;
U0S = 0.5;

flowParamsNames = {'BoxHeight','U0'};

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
nSamples=2000000;  

%nSamples=500000;  

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

Phys_Area = struct('shape','Box','N',[20;20],'L1',L1S,'L2',L2S);

Plot_Area = struct('y1Min',0,'y1Max',L1S,'N1',50,...
                       'y2Min',0,'y2Max',L2S,'N2',50);

V2_Num   = struct('Fex','Meanfield','N',[20;20],'L',1);

eq_Num   = struct('eqSolver','Newton','NewtonLambda1',0.1,'NewtonLambda2',1);
%eq_Num   = struct();

PhysArea = {Phys_Area};
PlotArea = {Plot_Area};


V2Num  =  {V2_Num};
eqNum  =  {eq_Num};

DDFTCode = {'DDFTDynamics'};
        
doPlots = false;

paramList = {'PhysArea','PlotArea','V2Num','eqNum','doPlots'};

DDFTParamsNames = {paramList, paramList, paramList, paramList};

% HIParamsNamesDDFT={'sigmaH','sigma'};               
HIParamsNamesDDFT={};
               
DDFTName={'r20'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r'};

% whether to do DDFT calculations
%doDDFT={true,true,true,true};
doDDFT={true};

% do we load and save the DDFT data
loadDDFT={true};

% for N computations
DDFTColour = {{'r','g','b'}};
DDFTStyle = {{'-','-','-'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

%plotType = 'contour';
plotType = 'surf';

% x axis for position and velocity plots
rMin=[0;0];
rMax=[L1S;L2S];

pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=2;

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
doMovieGif     = true;          % .gif movie
doInitialFinal = false;
doMeans        = true;
doEquilibria   = false;

%sendEmail = true;
