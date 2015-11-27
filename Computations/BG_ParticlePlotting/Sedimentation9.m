%--------------------------------------------------------------------------
% Physical setup
%--------------------------------------------------------------------------

% Names of potential files name.m to use for General, Initial, Final
% computations.  
% Potential .m files should go in the Potentials directory

% Note that the function should be of the form
% [V1,DV1]=potential(x,optsPhys) where
%   V1 is the one-body potential (vector of length(x)/dim), i.e. for
% each point of a DDFT vector, or each particle of a stochastic
% calculation.
%   DV1 is the gradient of V of length(x) used in stochastic calculation
V1DV1={'gravity','gravity','gravity'};
V2DV2='hardSphere';

% Names for the potentials to be used in creating the save files.
% Should be used to ensure we don't overwrite other data, and to make it
% easier to find the calculation later
% Of the form {General, Initial, Final}
potNames={'gravity','fixed9','gravity'};

% appropriate physical parameters for potentials in V1DV1
G={1,0,1};        % gravity strength
sep={0,2,0};    % separation of particles in initial state
% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParams=struct('G',G,'sep',sep);

% description for naming the save files
potDescr=cell(1,3);
for iDescr=1:3
    potDescr{iDescr}=[potNames{iDescr} '-' num2str(G{iDescr}) '-' num2str(sep{iDescr})];
end

% geometry; currently either 'planar' or 'spherical'.  
% Affects the way the 1D coordinates are calculated from the full 3D
% stochastic dynamics, and how the density is calculated for plots.
geom='planar';

% dimension in which to do the stochastic calculations
stocDim=3;
% dimension in which to do the DDFT calculations (=1)
DDFTDim=1;

% end time of calculation
tMax=40;

% physical constants:

m1=1;  % mass of first species
m2=2;  % mass of second species

sigma1=1;  % diameter of first species
sigma2=1;  % diameter of second species

sigmaH1=sigma1;  % hydrodynamic diameter of first species
sigmaH2=sigma2;    % hydrodynamic diameter of second species

nParticles1=5;  % number of first species
nParticles2=4;  % number of second species

nParticles=nParticles1+nParticles2;  % number of particles

m=[ m1*ones(nParticles1*stocDim,1) ; m2*ones(nParticles2*stocDim,1) ]; 
% column vector of masses by coordinate

sigma=[ sigma1*ones(nParticles1*stocDim,1) ; sigma2*ones(nParticles2*stocDim,1) ]; 
% column vector of diameters by coordinate

sigmaH=[ sigmaH1*ones(nParticles1*stocDim,1) ; sigmaH2*ones(nParticles2*stocDim,1) ]; 
% column vector of hydrodynamic diameters by coordinate

kBT=1;          % temperature

% Stokes-Einstein diffusion
eta=1/3/pi;
D0 = kBT./(3*pi*eta*sigmaH);  
% vector of diffusion coefficients by coordinate
gamma=kBT./m./D0;             
% vector of friction coefficients by coordinate

% could specify D0 and gamma in a different way if required


%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% no noise - means we only need to do one run
noise=false;

% do we have a fixed initial condition (true) or sample from the
% equilibrium distribution of the second potential in V1DV1 above
fixedInitial=true;

% function to determine:
% if fixedInitial=false: initial guess in slicesample
% if fixedInitial=true: the fixed initial condition
% should be of the form x = IC(optsPhys) and return a column vector of
% length nParticles*dim
initialGuess='NineParticles';

% do we sample the final equilibrium distribution or not.
% You don't want to do this unless you're evolving to an equilibrium, as
% otherwise its relationship to the calculation is pretty meaningless.
sampleFinal=false;

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true
nSamples=1;  
thin=1;             % see doc slicesample
burnin=0;        % see doc slicesample

% whether to load saved data for sampling
%loadSamples=false;
loadSamples=true;

% number of runs of stochastic dynamics to do, and average over
nRuns=1;

% the nRuns inital vectors are taken from the nSamples samples via
% runIC(i) = sample(sampleShift+sampleSkip*i)
% Note: these are related to burnin and thin, so you can either sample
% efficiently or choose efficiently from a more general sample.
% Defaults should be sampleShift=0 and sampleSkip=1
sampleShift=0;
sampleSkip=1;

% number of cores to use in parallel processing
poolsize=1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'rv','r','rv','r'};
% HI interaction matrices
stocHIType={[],[],'RPInv','RP'};

% whether to include hydrodynamic interactions
stocHI={false,false,true,true};

% names for stochastic calculations -- used as legend text
stocName={'Langevin, HI=0','EM, HI=0','Langevin, HI=1','EM, HI=1'};
% line styles for plotting stochastic data, note that the space are important
stocStyle={'  b','--b','  g','--g'};

% whether to do Langevin and Brownian dynamics
doStoc={true,false,false,false};
%doStoc={false,false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true,true};
%loadStoc={false,false,false,false};

% number of time steps
tSteps={10^4,10^4,10^4,10^4};

% whether to save output data (you probably should)
saveStoc=true;

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

% length of grid
L=8;
% number of points
%N=201;

N=201;
% parameter used in reallocating chebyshev points
ep=0.25;

% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'rv','r','rv','r'};

% whether to include hydrodynamic interactions
DDFTHI={false,false,true,true};

% names for DDFT calculations -- used as legend text
DDFTName={'DDFT r/v, HI=0','DDFT r, HI=0','DDFT r/v, HI=1','DDFT r, HI=1'};
% line styles for DDFT data, note that the space are important
DDFTStyle={'  r','--r','  m','--m'};

% whether to do DDFT calculations
%doDDFT={true,true,false,false};
doDDFT={false,false,false,false};

% do we load and save the DDFT data
loadDDFT={true,true,true,true};
%loadDDFT={false,false,false,false};

saveDDFT=true;

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

% plotting options in the form {General, Initial, Final}

% number of plots, which is also the number of times at which to save the
% stochastic run data
nPlots=100;
%nPlots=11;

% times at which to save and plot the stochastic simulations
plotTimes=0:tMax/nPlots:tMax;

% x axis for position and velocity plots
rMin={-1250,0,-1250};
rMax={10,10,-1235};
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin={0,0,0};
RMax={7,7,7};
PMin={-50,-10,-50};
PMax={0,0,-49};

% y axis for mean position and velocity plots
RMMin=-1250;
RMMax=10;
PMMin=-60;
PMMax=10;

% position of the legend
%legPos='SouthEast';
legPos='East';

% number of bins for histograming of stochastic data
nBins=50;

% don't plot v for regions where rho<vCutoff*max(rho) as these are just
% noise either due to DDFT calculation or very small number of particles
vCutoff=5*10^(-2);

% colours in which to plot the spheres
%colours={{'r','b'}};
colours={{['r','b','b','b','b','g','g','g','g']}};
% how much to shift each plot by (column vector for each calculation)
shift=[0; 0; 0];
% viewpoint of final 3D plot
%viewPoint=[-0.3 -1 -0.6];
viewPoint=[163 16];
% whether to renormalize each axis, i.e. whether to follow it in time
renormalize=[false; false; true];
% which particle to follow
%followType='max';
followType='com';
comPlot=1;
% whether to relabel axis ticks
relabel=[true; true; true];
% axis limits
lims={{[-4,4],[-4, 4],[-5,5]}};
% axis tick marks and corresponding labels
ticks={{[],[],[]}};
labs={{{''},{''},{''}}};
% colour for the bath (background)
bath=[0.9,1,1];
% scale for velocities so that arrows in movie are scaled independently of
% time
pScale=5;

% animation options:
% resolution
dpi=150;
% framerate
fps=5;
% whether to flatten pdfs in swf
bitmap=true;
% whether to suppress output of gs and pdf2swf
quiet=true;

% determine which movies/plots to make
doMovieGif=false;
doPdfs=false;
doMovieSwf=false;
doInitialFinal=false;
doMeans=false;

doInitialFinalP=false;
doMovieGifP=true;
doPdfsP=false;
doMovieSwfP=false;

doCustomP=false;
custom='ParticlesWithPath';
