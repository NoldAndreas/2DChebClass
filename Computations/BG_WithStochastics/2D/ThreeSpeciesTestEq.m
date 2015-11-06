%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

%nParticlesS=[20;20;20];
nParticlesS=[5;5;5];

kBT=1;          % temperature
mS=[1;1;1];

%mS = 1;
%gammaS  = 1;

gammaS=[1;1;1];

D0S=kBT./mS./gammaS;


%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='V1_Test_Box';

% appropriate physical parameters for potentials in V1DV1

y0 = 10;

% box
L1S = y0;
L2S = y0;
kBTS = kBT;

% potential
V0S = 2;
y10S = 7;
y20S = 6;
tauS = 10;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'V0','tau','y10','y20','L1','L2'};

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
tMax=7;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true
%nSamples=50000;  

nSamples=10000;  

sampleFinal = true;

initialGuess='makeGridPosScale';

% number of runs of stochastic dynamics to do, and average over
nRuns=100;

% number of cores to use in parallel processing
%poolsize=12;
poolsize=4;

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
%doStoc={false,false,false,false};
doStoc={true,false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={false,true,true,true};

% number of time steps
tSteps={10^3,10^3,2*10^4,10^3};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};

stocColour = {{'r','b','g'},{'g'},{'b'},{'m'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

Phys_Area = struct('shape','Box','N',[20,20],'y1Min',0,'y1Max',y0,...
                   'y2Min',0,'y2Max',y0);
Plot_Area = struct('y1Min',0,'y1Max',y0,'N1',100,...
                       'y2Min',0,'y2Max',y0,'N2',100);
V2Num   = struct('Fex','Meanfield','N',[20;20],'L',1);
               

PhysArea = {Phys_Area, Phys_Area};
PlotArea = {Plot_Area, Plot_Area};
V2Num   = {V2Num, V2Num};

DDFTCode = {'DDFTDynamics', ...
            'DDFTDynamics'};
        
doPlots = false;

DDFTParamsNames = {{'PhysArea','PlotArea','V2Num','doPlots'}, ...
                   {'PhysArea','PlotArea','V2Num','doPlots'}};

% HIParamsNamesDDFT={'sigmaH','sigma'};               
HIParamsNamesDDFT={};
               
DDFTName={'r0','r1'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r','r'};

% whether to do DDFT calculations
doDDFT={false,false};
%doDDFT={false,false};

% do we load and save the DDFT data
loadDDFT={true,true};

DDFTColour = {{'g','m','c'},{'r','b'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

% x axis for position and velocity plots
rMin=[-2;-2];
rMax=[y0+2;y0+2];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.2;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[0;0];
RMMax=[y0;y0];
PMMin=[-1;-1];
PMMax=[1;1];

% number of bins for histograming of stochastic data
nBins=[50;50];

% determine which movies/plots to make
% distribution movies/plots
doMovieGif=false;          % .gif movie
doInitialFinal=false;
doMeans=false;
doEquilibria = true;

%sendEmail = true;
