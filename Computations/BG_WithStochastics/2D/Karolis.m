%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

%nParticlesS=50;
nParticlesS=60;

kBT=1;          % temperature
mS=1;

gammaS=1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='free2D_Box';

% appropriate physical parameters for potentials in V1DV1
L1S = 10;
L2S = 10;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'L1','L2'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='hardSphere';

sigmaS = 1;

potParams2Names={'sigma'};

%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------

sigmaHS = 1;

HIParamsNames={'sigmaH'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=0.1;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true

nSamples=2000000;  
burnin = 50000;

initialGuess='makeGridPosScale';

% number of runs of stochastic dynamics to do, and average over

nRuns=1;

% number of cores to use in parallel processing
%poolsize=12;
poolsize=1;

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
doStoc={true,false,true,false};
%doStoc={false,false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true,true};

% number of time steps
tSteps={10^4,10^3,10^4,10^3};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};

stocColour = {{'r'},{'g'},{'b'},{'m'}};



%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

viewPoint = [-56;7];

% x axis for position and velocity plots
rMin=[0;0];
rMax=[L1S;L2S];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=1;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[-3;5];
RMMax=[3;7];
PMMin=[-1;-1];
PMMax=[1;1];

% number of bins for histograming of stochastic data
nBins=[40;40];

% determine which movies/plots to make
% distribution movies/plots
doMovieGif     = false;          % .gif movie
doInitialFinal = false;
doMeans        = false;
doEquilibria   = true;

%sendEmail = true;
