%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

nParticlesS=10;

kBT=1;          % temperature
mS=1;

gammaS=1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='infHIDiffusion3';

% appropriate physical parameters for potentials in V1DV1
V0S = 0.01;
V0addS = 3;
tauS = 0.1;
sigma1AddS = 0.5;
sigma2AddS = 0.5;
y10S = -1;
y20S = -1;
y11S = 1;
y21S = -1;
y12S = 0;
y22S = 0.5;



% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'V0','V0add','tau','sigma1Add','sigma2Add',...
                  'y10','y20','y11','y21','y12','y22'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

% V2DV2='hardSphere';
% 
% sigmaS = 1;
% 
% potParams2Names={'sigma'};

V2DV2='Gaussian';

epsilonS = 1;
alphaS   = 1;

potParams2Names={'epsilon','alpha'};



%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------

sigmaHS = 0.5;

HIParamsNames={'sigmaH'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=0.25;

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
nRuns=15000;

% number of cores to use in parallel processing
poolsize=12;
%poolsize=1;

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
loadStoc={true,true,true,true};

% number of time steps
tSteps={10^3,10^3,2*10^4,10^3};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};

stocColour = {{'r'},{'g'},{'b'},{'m'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

y0 = 3;

PhysArea = {struct('y1Min',-inf,'y1Max',inf,'N',[30,30],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4), ...
            struct('y1Min',-inf,'y1Max',inf,'N',[30,30],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4)};
PlotArea = {struct('y1Min',-y0,'y1Max',y0,'N1',100,...
                       'y2Min',-y0,'y2Max',y0,'N2',100), ...
            struct('y1Min',-y0,'y1Max',y0,'N1',100,...
                       'y2Min',-y0,'y2Max',y0,'N2',100)};
% FexNum   = {struct('Fex','FMTRosenfeld',...
%                        'Ncircle',10,'N1disc',10,'N2disc',10), ...
%             struct('Fex','FMTRosenfeld',...
%                        'Ncircle',10,'N1disc',10,'N2disc',10)};

FexNum  = {struct('Fex','Meanfield','N',[20;20],'L',2), ...
            struct('Fex','Meanfield','N',[20;20],'L',2)};

 
% HINum    = {[], ...
%             struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','RP12_2D', ...
%                       'HIPreprocess', 'RotnePragerPreprocess2D', ...
%                       'sigma',sigmaS,'sigmaH',sigmaS/2)};
                  

HINum = {[],[]};

DDFTCode = {'DDFT_DiffusionInfSpace_NSpecies2', ...
            'DDFT_DiffusionInfSpace_NSpecies2'};
        
Tmax = tMax;

doPlots = true;

DDFTParamsNames = {{'PhysArea','PlotArea','FexNum','Tmax','doPlots'}, ...
                   {'PhysArea','PlotArea','FexNum','HINum','Tmax','doPlots'}};

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

DDFTColour = {{'g'},{'m'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

% whether to plot the distribution (false) or the density (true) in
% spherical coordinates.
plotDensity=false;

plotCurrent=false;

plotType = 'surf';

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

viewPoint = [-30,45];

% position of the legend -- set to 'off' to turn off individual legends
%legPos='SouthEast';
%legPos='NorthWest';
legPos='off';
% to create a single legend for all plots
oneLeg='top';
%oneLeg='off';
perRow=4;

% number of bins for histograming of stochastic data
nBins=[20;20];
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
doMeans=false;

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
