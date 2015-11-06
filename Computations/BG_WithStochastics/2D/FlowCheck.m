%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

nParticlesS=2;

kBT=1;          % temperature
mS=1;

gammaS=1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='free2D';

potParamsNames = {};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='Gaussian';

epsilonS = 0;
alphaS   = 1;

potParams2Names={'epsilon','alpha'};

%--------------------------------------------------------------------------
% Flow parameters
%--------------------------------------------------------------------------

UDU = 'FlowCheckU';

U0S = 1;
flowParamsNames = {'U0'};

%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------

sigmaHS = 0.5;

HIParamsNames={'sigmaH'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=30;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

fixedInitial = true;

noise = false;

initialGuess = [2,0,1,1];

% number of runs of stochastic dynamics to do, and average over
nRuns=1;
nSamples = 1;

% number of cores to use in parallel processing
poolsize=1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r'};

% whether to include hydrodynamic interactions
stocHI={false};
% HI interaction matrices
stocHIType={[]};

% names for stochastic calculations -- used as legend text
stocName={'testFlow'};

% whether to do Langevin and Brownian dynamics
doStoc={true};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={false};

% number of time steps
tSteps={10^4};

% whether to save output data (you probably should)
saveStoc={true};

stocColour = {{'r'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

y0 = 5;

% x axis for position and velocity plots
rMin=[-y0;-y0];
rMax=[y0;y0];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.2;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[-y0;-y0];
RMMax=[y0;y0];
PMMin=[-1;-1];
PMMax=[1;1];

% number of bins for histograming of stochastic data
nBins=[30;30];

% determine which movies/plots to make
% distribution movies/plots
doMovieGif=false;          % .gif movie
doInitialFinal=false;
doMeans=false;



%sendEmail = true;
