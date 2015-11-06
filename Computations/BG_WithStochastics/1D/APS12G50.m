%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='spherical';

% dimension in which to do the stochastic calculations
stocDim=3;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=1;

nParticlesS=[25;25];
%nParticlesS=[10;10];

kBT=1;          % temperature
mS=[1;1];

gammaS=[2;2];
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='APSG';

% appropriate physical parameters for potentials in V1DV1   
VmS= [0.01;0.01];  
tSwitchS= [12;12];

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'Vm','tSwitch'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='Gaussian';

epsilonS = [0.5,0.5;0.5,0.5];
%epsilonS = [0,0 ; 0,0];
alpha1=0.2;
alpha2=1;
alpha12=(alpha1+alpha2)/2;

alphaS = [alpha1, alpha12; alpha12, alpha2];

potParams2Names={'epsilon','alpha'};

%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------
% 
% HIParamsNames={};
% HIParamsNamesDDFT={};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=12;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true
%nSamples=50000;  

nSamples=50000;  

initialGuess='makeGrid';

% number of runs of stochastic dynamics to do, and average over
nRuns=1000;

% number of cores to use in parallel processing
poolsize=12;
%poolsize=1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r','rv'};

% whether to include hydrodynamic interactions
stocHI={false,false};
% HI interaction matrices
stocHIType={[],[]};

% names for stochastic calculations -- used as legend text
stocName={'r0','rv0'};

% whether to do Langevin and Brownian dynamics
doStoc={true,true};
%doStoc={false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true};
%loadStoc={false,false};

% number of time steps
tSteps={10^3,10^3};

% whether to save output data (you probably should)
saveStoc={true,true};


%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

Phys_Area = struct('shape','InfSpectralLineSpherical','N',200,'L',4);
Plot_Area = struct('N',200,'yMin',0,'yMax',10);
Fex_Num   = struct('Fex','Meanfield','N',100,'L',2);

PhysArea = {Phys_Area,Phys_Area};
PlotArea = {Plot_Area,Plot_Area};
FexNum   = {Fex_Num,Fex_Num};

DDFTCode = {'DDFT_Diffusion_1D_Spherical', ...
            'DDFT_Inertia_1D_Spherical'};
        
doPlots = true;

DDFTParamsNames = {{'PhysArea','PlotArea','FexNum','doPlots'}, ...
                   {'PhysArea','PlotArea','FexNum','doPlots'}};
                           
DDFTName={'r0','rv0'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r','rv'};

% whether to do DDFT calculations
doDDFT={true,true};

% do we load and save the DDFT data
loadDDFT={true,true};


%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------


% whether to plot the distribution (false) or the density (true) in
% spherical coordinates.
plotDensity=false;

plotCurrent=false;

% x axis for position and velocity plots
rMin=0;
rMax=12;
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=15;

PMin=-0.5;
PMax=0.4;


% y axis for mean position and velocity plots
RMMin=1;
RMMax=4;
PMMin=-0.5;
PMMax=0.25;

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
doInitialFinal=true;
doMeans=true;

% particle movies/plots
doInitialFinalP=false;
doMovieGifP=false;
doPdfsP=false;
doMovieSwfP=false;

% strip out any where the time step was too large
doStrip=false;

doCustom=false;
%custom='makeMovieWithSnapshotsV';
custom='MSHISnapshots';
