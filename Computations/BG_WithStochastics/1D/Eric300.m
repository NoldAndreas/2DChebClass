%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='spherical';

% dimension in which to do the stochastic calculations
stocDim=3;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=1;

nParticlesS=300;

kBT=1;          % temperature
mS=1;

r0 = 1;%0.080961;

%gammaS=6*pi*r0;
gammaS = 1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% IC parameters
%--------------------------------------------------------------------------

ICrho = 'EricIC';

ICParamsNames = {'nParticles'};

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='EricV1';

% appropriate physical parameters for potentials in V1DV1
VmS = 0.301;    
VqS = 0.328;
y0S = 31;
y1S = 25;

alphaS = 2.47;
%alphaS = 1;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'Vm','Vq','y0','y1','alpha'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='EricV2';

sigmaS = r0;  % particle radius
%V0S    = 1;%100;

potParams2Names={'sigma'};%,'V0'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=1;

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------


Phys_Area = struct('shape','InfSpectralLineSpherical','N',200,'L',20);
Plot_Area = struct('N',200,'yMin',0,'yMax',40);
Fex_Num   = struct('Fex','FMT','N',100);

PhysArea = {Phys_Area};
PlotArea = {Plot_Area};
FexNum   = {Fex_Num};
HINum    = {[]};

DDFTCode = {'DDFT_Inertia_1D_Spherical'};
%DDFTCode = {'DDFT_Diffusion_1D_Spherical_GivenIC'};
       
doPlots = true;

DDFTParamsNames = {{'PhysArea','PlotArea','FexNum','doPlots'}};
               
%DDFTName={'Inertial DDFT'};
DDFTName={'Diffusive DDFT'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
%DDFTType={'rv'};
DDFTType={'r'};

% whether to do DDFT calculations
doDDFT = {true};

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
rMin=0;
rMax=5;
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=300;

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

% determine which movies/plots to make
% distribution movies/plots
doMovieGif=true;          % .gif movie
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
%custom='makeMovieWithSnapshotsV';
custom='MSHISnapshots';
