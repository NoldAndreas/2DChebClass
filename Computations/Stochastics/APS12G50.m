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

kBT=1;          % temperature
mS=[1;1];

gammaS=[2;2];
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1={'APSG'};

% appropriate physical parameters for potentials in V1DV1   
VmS= 0.01;  

tSwitchS=12;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {{'Vm','tSwitch'}};

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

potParams2Names={{'epsilon','alpha'}};

%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------

HIParamsNames={{}};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=12;

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

% length of grid
L=4;
% number of points
N=201;

% plotting range of the form (rPlotMin:rPlotStep:rPlotMax);
rPlotMin =0;
rPlotMax=10;
rPlotStep=0.01;

intLimMin=0;
intLimMax=20;

NInt=N;

% integration region
uBound=2;
lBound=0;


FexMatrices='getRPAMatricesFull';
Fex='RPA';

% file names for the DDFT integration matrices
DDFTHIType={'noHI','noHI','JeffreyOnishi11','RotnePrager'};

DDFTParamsNames={{'L','N','rPlotMin','rPlotMax','rPlotStep', ...
                   'intLimMin','intLimMax','NInt','uBound','lBound', ...
                   'FexMatrices','Fex','DDFTHIType'}};

DDFTParamsSaveNames={{'N'}};

DDFTName={'DDFT 1','DDFT 2','DDFT 3','DDFT 4'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'rv','r','rv','r'};

% whether to include hydrodynamic interactions
% currently this should be false
DDFTHI={false,false,true,true};

% locations for DDFT code
DDFTDir=['DDFT1D' filesep 'Spherical'];
DDFTCode='DDFTSpherical';

% whether to do DDFT calculations
%doDDFT={true,true,true,true};
doDDFT={true,false,false,false};
%doDDFT={false,false,false,false};

% do we load and save the DDFT data
%loadDDFT={true,true,true,true};
loadDDFT={false,false,false,false};

saveDDFT=true;


%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

% whether to plot the distribution (false) or the density (true) in
% spherical coordinates.
plotDensity=false;

plotCurrent=false;

% x axis for position and velocity plots
% in form {General (movie), Initial, Final}
rMin={0,0,0};
rMax={10,10,10};
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
% in form {General (movie), Initial, Final}
RMin={0,0,0};

if(plotDensity)
    RMax={1.5,1.5,1.5};        % density
else
    %RMax={100,100,100};    % distribution
    RMax={15,15,15};    % distribution
end

PMin={-0.6,-0.6,-0.6};  % velocity
PMax={0.4,0.4,0.4};
% PMin={-20,-20,-20};  % current
% PMax={10,10,10};


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
