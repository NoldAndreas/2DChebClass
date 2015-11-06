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

V1DV1='APSHS';

% appropriate physical parameters for potentials in V1DV1
VmS= 0.1;    

tSwitchS=1;

ZS=[0.5;0.25];

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'Vm','Z','tSwitch'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='hardSphere';

sigma1=1;
sigma2=1;
sigma12=(sigma1+sigma2)/2;

sigmaS=[sigma1  sigma12;
        sigma12 sigma2];

potParams2Names={'sigma'};

%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------

sigmaH1=0.5;
sigmaH2=0.5;
sigmaH12=(sigmaH1+sigmaH2)/2;

sigmaHS=[sigmaH1  sigmaH12;
        sigmaH12 sigmaH2];

HIParamsNames={'sigmaH'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=2.5;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true
nSamples=50000;  

initialGuess='makeGrid';

% number of runs of stochastic dynamics to do, and average over
%nRuns=2;

nRuns = 10;

% number of cores to use in parallel processing
%poolsize=12;
poolsize=4;

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
%doStoc={true,true,false,false};
doStoc={true,false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true,true};
%loadStoc={false,false,false,false};

% number of time steps
tSteps={2*10^4,2*10^4,2*10^4,2*10^4};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};


%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------


Phys_Area = struct('shape','InfSpectralLineSpherical','N',200,'L',4);
Plot_Area = struct('N',200,'yMin',0,'yMax',12);
Fex_Num   = struct('Fex','FMT','N',100);

PhysArea = {Phys_Area, Phys_Area, Phys_Area, Phys_Area, Phys_Area};
PlotArea = {Plot_Area, Plot_Area, Plot_Area, Plot_Area, Plot_Area};
FexNum   = {Fex_Num, Fex_Num, Fex_Num, Fex_Num, Fex_Num};
HINum    = {[], ...
            [], ...
            struct('N',100,'L',4, ...
                       'HI11','noHISpherical',...
                       'HI12','RotnePrager12Spherical', ...
                       'HIPreprocess', 'RotnePragerPreprocessSpherical'), ...
            struct('N',100,'L',4, ...
                       'HI11','JeffreyOnishi11Spherical',...
                       'HI12','JeffreyOnishi12Spherical', ...
                       'HIPreprocess', 'JeffreyOnishiPreprocessSpherical'), ...
            struct('N',100,'L',4, ...
                       'HI11','JeffreyOnishi11Spherical',...
                       'HI12','JeffreyOnishi12Spherical', ...
                       'HIPreprocess', 'JeffreyOnishiPreprocessSpherical', ...
                       'HIg','g_0_Spherical')};

DDFTCode = {'DDFT_Diffusion_1D_Spherical', ...
            'DDFT_Inertia_1D_Spherical', ...
            'DDFT_Diffusion_1D_Spherical', ...
            'DDFT_Inertia_1D_Spherical', ...
            'DDFT_Inertia_1D_Spherical'};
       
doPlots = false;

DDFTParamsNames = {{'PhysArea','PlotArea','FexNum','doPlots'}, ...
                   {'PhysArea','PlotArea','FexNum','doPlots'},...
                   {'PhysArea','PlotArea','FexNum','HINum','doPlots'},...
                   {'PhysArea','PlotArea','FexNum','HINum','doPlots'}, ...
                   {'PhysArea','PlotArea','FexNum','HINum','doPlots'}};

HIParamsNamesDDFT={'sigmaH','sigma'};               
               
DDFTName={'r0','rv0','r1','rv1','rv1gTest'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r','rv','r','rv','rv'};

% whether to do DDFT calculations
%doDDFT={true,true,true,true,true};
doDDFT={true,false,false,false,false};
%doDDFT={false,false,false,false,false};

% do we load and save the DDFT data
loadDDFT={true,true,true,true,true};
%loadDDFT={false,false,false,false,false};

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

% determine which movies/plots to make
% distribution movies/plots
doMovieGif=false;          % .gif movie
doPdfs=false;              % .pdfs to make .swf
doMovieSwf=false;          % .swf movie
doInitialFinal=false;
doMeans=false;
doEquilibria = true;

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
