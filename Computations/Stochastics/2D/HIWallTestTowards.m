%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

nParticlesS=20;

kBT=1;          % temperature
mS=1;

gammaS=1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='V1_Well_gravity';

% appropriate physical parameters for potentials in V1DV1
V0S        = 0.01;
V0addS     = 3;
tauS       = 0.1;
%sigma1AddS = 0.5;
%sigma2AddS = 0.5;

sigma1AddS = 1;
sigma2AddS = 1;

y10S       = 0;
y20S       = 2;
gS         = 1;
gcutS      = 10; 
% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'V0','V0add','tau','sigma1Add','sigma2Add',...
                  'y10','y20','g','gcut'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='hardSphere';

sigmaS = 1;

potParams2Names={'sigma'};

%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------

sigmaHS = 0.5;

HIParamsNames={'sigmaH'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=0.2;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true

nSamples=500000;  
 

initialGuess='makeGrid';

% number of runs of stochastic dynamics to do, and average over

nRuns=10000;

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
%doStoc={true,false,true,false};
doStoc={false,false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true,true};

% number of time steps
tSteps={10^4,10^3,10^4,10^3};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};

stocColour = {{'r'},{'g'},{'b'},{'m'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

Phys_Area = struct('shape','HalfSpace_FMT','N',[25;25],'L1',2,'L2',2, ...
                       'y2wall',0,'N2bound',10,'h',1,'L2_AD',1,'alpha_deg',90); 

Sub_Area = struct('shape','Box','y1Min',-3,'y1Max',3,'N',[20,20],...
                      'y2Min',0.5,'y2Max',1);
                   
Plot_Area = struct('y1Min',-3,'y1Max',3,'N1',100,...
                       'y2Min',0.5,'y2Max',5,'N2',100);

Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);
                   
PhysArea = {Phys_Area, Phys_Area, Phys_Area, Phys_Area};

SubArea  = {Sub_Area, Sub_Area, Sub_Area, Sub_Area};

PlotArea = {Plot_Area, Plot_Area, Plot_Area, Plot_Area};

FexNum   = {Fex_Num, Fex_Num, Fex_Num, Fex_Num};

V2Num    = {[],[],[],[]};

HINum    = {[], ...
            struct('N',[30;30],'L',2,'HI11','noHI_2D','HI12','FullWallHI_2D_noConv', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D',...
                      'HIWallFull',true,'doConv',false,...
                      'Wall','SelfWallTermJK'), ...
            struct('N',[30;30],'L',2,'HI11','noHI_2D','HI12','noHI_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D',...
                      'HIWallFull',true,'doConv',false,...
                      'Wall','SelfWallTermJK'), ...
            struct('N',[30;30],'L',2,'HI11','noHI_2D','HI12','Oseen_2D_noConv', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D',...
                      'HIWallFull',true,'doConv',false,...
                      'Wall','SelfWallTermJK'), ...
           };

DDFTCode = {'DDFTDynamics', 'DDFTDynamics', 'DDFTDynamics', 'DDFTDynamics'};
        
doPlots = true;

DDFTParamsNames = {{'PhysArea','SubArea','PlotArea','FexNum','V2Num','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','HINum','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','HINum','doPlots'}};

HIParamsNamesDDFT={'sigmaH','sigma'};               
               
DDFTName={'No HI','Full HI','Just Wall','Oseen + Wall'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r','r','r','r'};

% whether to do DDFT calculations
doDDFT={true,true,false,false};

% do we load and save the DDFT data
loadDDFT={true,true,true,true};
%loadDDFT={false,false};

DDFTColour = {{'r'},{'b'},{'g'},{'m'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

viewPoint = [-56;7];

% x axis for position and velocity plots
rMin=[-3;0];
rMax=[3;5];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.4;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[-3;5];
RMMax=[3;7];
PMMin=[-1;-1];
PMMax=[1;1];

% number of bins for histograming of stochastic data
nBins=[20;20];

% determine which movies/plots to make
% distribution movies/plots
doMovieGif     = true;          % .gif movie
doInitialFinal = true;
doMeans        = false;
doEquilibria   = false;

%sendEmail = true;
