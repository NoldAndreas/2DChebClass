%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='full';

% dimension in which to do the stochastic calculations
stocDim=3;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=1;

nParticlesS=9;

kBT=1;          % temperature
mS=1;

gammaS=1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='gravity3D';

% appropriate physical parameters for potentials in V1DV1
gS= 1;    

sepS = 2;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'g','sep'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='hardSphere';

sigmaS=1;

potParams2Names={'sigma'};

%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------

sigmaHS = sigmaS;

HIParamsNames={'sigmaH'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=40;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

noise = true;
fixedInitial = true;

initialGuess='NineParticles';

sampleFinal = false;

nSamples = 1;
nRuns = 1;

poolsize = 1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r','rv','r','rv'};

% use new hard sphere code with resampling for collisions
useNewHS = {true,false,false,false};

% whether to include hydrodynamic interactions
stocHI={false,false,true,true};
% HI interaction matrices
stocHIType={[],[],'RP','RPInv'};

% names for stochastic calculations -- used as legend text
stocName={'r0','rv0','r1','rv1'};

% whether to do Langevin and Brownian dynamics
%doStoc={true,true,true,true};
%doStoc={true,true,false,false};
doStoc={true,false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={false,true,true,true};
%loadStoc={false,false,false,false};

% number of time steps
tSteps={10^4,10^4,10^4,10^4};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};


%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

%nPlots = 10;
%plotTimes=0 : tMax / nPlots : tMax;

% x axis for position and velocity plots
rMin=0;
rMax=12;
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=12;

PMin=-0.5;
PMax=0.2;


% y axis for mean position and velocity plots
RMMin=3.4;
RMMax=4.1;
PMMin=-0.3;
PMMax=0.1;

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


% colours in which to plot the spheres
%colours={{'r','b'}};
colours={{['r','b','b','b','b','g','g','g','g'],['r','b','b','b','b','g','g','g','g']}};
% how much to shift each plot by (column vector for each calculation)
%shift=[0, 5; 0, 0; 0, 0];
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
%lims={{[-10,10],[-4, 4],[-5,5]}};
lims={{[-10,10],[-10, 10],[-10,10]}};
% axis tick marks and corresponding labels
ticks={{[],[],[]}};
labs={{{''},{''},{''}}};
% colour for the bath (background)
bath=[0.9,1,1];
% scale for velocities so that arrows in movie are scaled independently of
% time
pScale=5;




% particle movies/plots
doInitialFinalP=false;
doMovieGifP=false;
doPdfsP=true;
doMovieSwfP=false;
