%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar';

% dimension in which to do the stochastic calculations
stocDim=1;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=1;

nParticlesS=10;

kBT=1;          % temperature
mS=1;

%gammaS=0.1; % fails
gammaS=0.2; % fails
%gammaS=0.3; % fails
%gammaS=0.5; % fails
%gammaS=1; % fails
%gammaS=2; % runs but messy

D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='quadraticGaussianFade21D';

% appropriate physical parameters for potentials in V1DV1
VmS        = 0.01;    
V0AddS     = 15;
sigma2AddS = 10;
y0S        = 4;
tauS       = 0.01;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'Vm','V0Add','sigma2Add','y0','tau'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='Gaussian';

epsilonS = 0;  % no interaction
alphaS   = 1;

potParams2Names={'epsilon','alpha'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=5;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true
nSamples = 50000;
burnin = 10000;

initialGuess='makeGrid';

% number of runs of stochastic dynamics to do, and average over
%nRuns=2;

nRuns = 50000;

%nRuns = 50000;

% number of cores to use in parallel processing
poolsize=12;
%poolsize=4;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'rv'};

% whether to include hydrodynamic interactions
stocHI={false};
% HI interaction matrices
stocHIType={[]};

% names for stochastic calculations -- used as legend text
stocName={'Stoc'};

% whether to do Langevin and Brownian dynamics
doStoc={true};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true};

% number of time steps
tSteps={10^3};

% whether to save output data (you probably should)
saveStoc={true};


%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------


Phys_Area = struct('shape','InfSpectralLine','N',200,'L',4);
Plot_Area = struct('N',200,'yMin',-10,'yMax',10);
Fex_Num   = struct('Fex','Meanfield','N',100,'L',2);

PhysArea = {Phys_Area};
PlotArea = {Plot_Area};
FexNum   = {Fex_Num};

DDFTCode = {'DDFT_Inertia_1D_Planar'};
        
doPlots = true;

DDFTParamsNames = {{'PhysArea','PlotArea','FexNum','doPlots'}};
                           
DDFTName={'DDFT'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'rv'};

% whether to do DDFT calculations
doDDFT={false};

% do we load and save the DDFT data
loadDDFT={true};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

% whether to plot the distribution (false) or the density (true) in
% spherical coordinates.
plotDensity=false;

plotCurrent=false;

% x axis for position and velocity plots
rMin=-10;
rMax=10;
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=12;

PMin=-0.5;
PMax=0.5;


% y axis for mean position and velocity plots
RMMin=-10;
RMMax=10;
PMMin=-1;
PMMax=1;

% position of the legend -- set to 'off' to turn off individual legends
%legPos='SouthEast';
%legPos='NorthWest';
legPos='off';
% to create a single legend for all plots
oneLeg='top';
%oneLeg='off';
perRow=4;

% number of bins for histograming of stochastic data
nBins=50;
%nBins=100;

% don't plot v for regions where rho<vCutoff*max(rho) as these are just
% noise either due to DDFT calculation or very small number of particles
%vCutoff=5*10^(-2);
vCutoff=5*10^(-2);

% animation options:
% resolution
dpi=300;
% framerate
fps=5;
% whether to flatten pdfs in swf
bitmap=false;
% whether to suppress output of pdf2swf
quiet=true;

% determine which movies/plots to make
% distribution movies/plots
doMovieGif=false;          % .gif movie
doPdfs=false;              % .pdfs to make .swf
doMovieSwf=false;          % .swf movie
doInitialFinal=false;
doMeans=false;
doEquilibria = false;

% particle movies/plots
doInitialFinalP=false;
doMovieGifP=false;
doPdfsP=false;
doMovieSwfP=false;

% strip out any where the time step was too large
doStrip=false;

doCustom=false;
