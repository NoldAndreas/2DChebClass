%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

nParticlesS=[5;5];

kBT=1;          % temperature
mS=[1;1];

gammaS=[1;1];
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

%V1DV1='free2D_Box';
V1DV1 = 'quad_Box';

% appropriate physical parameters for potentials in V1DV1
L1S = 8;
L2S = 6;
kBTS = kBT;

tauS = 1;
y10S = 5;
y20S = 5;
B10S = 1;
B20S = 1;
V0AddS  = 1;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'L1','L2', ...
                  'tau','y10','y20','B10','B20','V0Add'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='Gaussian';

epsilonS = [1 1; 1 1];
alphaS   = [1 1; 1 1];

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

nSamples=5000;  

initialGuess='makeGridPos';

% number of runs of stochastic dynamics to do, and average over
%nRuns=500000;
nRuns=10;

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
doStoc={true,false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true,true};

% number of time steps
tSteps={10^4,10^3,2*10^4,10^3};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};

stocColour = {{'r','b'},{'g'},{'b'},{'m'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

y0 = 3;

Phys_Area = struct('shape','Box','N',[30,30],'L1',L1S,'L2',L2S);
Plot_Area = struct('y1Min',0,'y1Max',L1S,'N1',100,...
                       'y2Min',0,'y2Max',L2S,'N2',100);
V2_Num   = struct('Fex','Meanfield','N',[20;20],'L',2);

PhysArea = {Phys_Area, Phys_Area};
PlotArea = {Plot_Area, Plot_Area};


V2Num  =  {V2_Num, V2_Num};

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
doDDFT={true,false};
%doDDFT={false,false};

% do we load and save the DDFT data
loadDDFT={true,true};

DDFTColour = {{'g','m'},{'m'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

% x axis for position and velocity plots
% rMin=[0;0];
% rMax=[L1S;L2S];
rMin=[-2;-2];
rMax=[L1S+2;L2S+2];


pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.5;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[0;0];
RMMax=[max(L1S,L2S);max(L1S,L2S)];
PMMin=[-1;-1];
PMMax=[1;1];

% number of bins for histograming of stochastic data
nBins=[30;30];

% determine which movies/plots to make
% distribution movies/plots
doMovieGif     = false;          % .gif movie
doInitialFinal = false;
doMeans        = false;
doEquilibria   = true;

%sendEmail = true;
