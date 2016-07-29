%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='full';

stocDim=3;
DDFTDim=3;

nParticlesS=100;

kBT=1;          % temperature
mS=1;

gammaS=1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='freeBox3D';

% appropriate physical parameters for potentials in V1DV1

% box
L1S = 10;
L2S = 10;
L3S = 100;
%kBTS = kBT;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'L1','L2','L3'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='hardSphere';

sigmaS = 1;

potParams2Names={'sigma'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax = 0;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true
nSamples=1000000;  % done
burnin = 100000;

initialGuess='makeGridPosScale3D';

% number of runs of stochastic dynamics to do, and average over
nRuns=0;

% number of cores to use in parallel processing
poolsize=1;
%poolsize=4;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r'};

% whether to include hydrodynamic interactions
stocHI={false};
% HI interaction matrices
stocHIType={[]};

% names for stochastic calculations -- used as legend text
stocName={'Sampling'};

% whether to do Langevin and Brownian dynamics
doStoc={true};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true};

% number of time steps
tSteps={1};

% whether to save output data (you probably should)
saveStoc={true};

stocColour = {{'r'}};
stocStyle = {{'-'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

%plotType = 'contour';
%plotType = 'surf';

% x axis for position and velocity plots

%for stochastic comparison
%rMin=[0;0];
%rMax=[L1S;L2S];

% for N computation
%rMin = [0;0];
%rMax = [10;10];


rMin = [4;3];
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
nBins=[30;30];

%viewPoint=[-37.5,30];

separateSpecies = true;
separateError = true;

% determine which movies/plots to make
% distribution movies/plots
doMovieGif     = false;          % .gif movie
doInitialFinal = false;
doMeans        = false;
doEquilibria   = false;

sendEmail = false;
