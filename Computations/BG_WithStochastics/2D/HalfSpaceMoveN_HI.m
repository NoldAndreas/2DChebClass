%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

%nParticlesS=10;
nParticlesS=50;

kBT=1;          % temperature
mS=1;

gammaS=1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='V1_Well_Move_HalfSpace';

% appropriate physical parameters for potentials in V1DV1
V0S        = 0.01;

nAddS = 2;

V0addS     = 3;
tauS       = 0.01;

sigma1AddS = 250;
sigma2AddS = 50;

y10aS       = 0;
y10bS       = 0;
y20aS       = 8;
y20bS       = 4;

% y20aS       = 10;
% y20bS       = 2;


% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'V0','V0add','tau','sigma1Add','sigma2Add',...
                  'y10a','y20a','y10b','y20b','nAdd'};

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
%tMax=1.5;
tMax=1;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true


nSamples = 1000000;  

initialGuess='makeGridPos';

sampleFinal = false;

% number of runs of stochastic dynamics to do, and average over


%nRuns = 5000;

nRuns = 10000;

%nRuns = 100;

% number of cores to use in parallel processing
poolsize=12;
%poolsize=1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r','r','r','r'};

% whether to include hydrodynamic interactions
stocHI={false,true,true,true};
stocUseDivergence = {false,true,true,true};


% HI interaction matrices
stocHIType={[],'fullWall2D','wallMobility2D','RP2D'};

% names for stochastic calculations -- used as legend text
stocName={'No HI','Full HI','Only Wall', 'RP'};

% whether to do Langevin and Brownian dynamics
doStoc={true,true,false,true};
%doStoc={false,false,false,true};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true,true};

% number of time steps
tSteps={5*10^4,5*10^4,5*10^4,5*10^4};

% whether to save output data (you probably should)
saveStoc={true,true,true,true};

stocStyle = {{'-'},{'-'},{'-'},{'-'}};
stocColour = {{'r'},{'g'},{'b'},{'m'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

y1Plot=10;
y2Plot=10;

Phys_Area = struct('shape','HalfSpace_FMT','N',[40;40],'L1',3,'L2',3, ...
                       'y2wall',0,'N2bound',10,'h',1,'L2_AD',1,'alpha_deg',90); 

% Phys_Area = struct('shape','HalfSpace_FMT','N',[20;20],'L1',3,'L2',3, ...
%                        'y2wall',0,'N2bound',10,'h',1,'L2_AD',1,'alpha_deg',90);

Sub_Area = struct('shape','Box','y1Min',-3,'y1Max',3,'N',[20,20],...
                      'y2Min',0.5,'y2Max',1);
                   
Plot_Area = struct('y1Min',-10,'y1Max',10,'N1',20,...
                       'y2Min',0.5,'y2Max',10,'N2',20);

Fex_NumRosenfeld   = struct('Fex','FMTRosenfeld',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

Fex_NumRoth   = struct('Fex','FMTRoth',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

Fex_Num3D   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

HI_None = [];
HI_Full = struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','FullWallHI_RP_2D_noConv', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D',...
                      'HIWallFull',true,'doConv',false,...
                      'Wall','SelfWallTermKN');
                  
HI_OnlyWall = struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','noHI_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D',...
                      'HIWallFull',true,'doConv',false,...
                      'Wall','SelfWallTermKN');

HI_RP = struct('N',[20;20],'L',4,'HI11','noHI_2D','HI12','RP12_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D');                 
                  
HINum    = {HI_None, ...
            HI_Full, ...
            HI_OnlyWall, ...
            HI_RP, ...
           };

HIParamsNamesDDFT={'sigmaH','sigma'};                  
                  
%eq_Num    = struct('eqSolver','Newton','NewtonLambda1',0.7,'NewtonLambda2',0.7);
eq_Num = struct('eqSolver','fsolve');
                   
PhysArea = {Phys_Area, Phys_Area, Phys_Area, Phys_Area};

SubArea  = {Sub_Area, Sub_Area, Sub_Area, Sub_Area};

PlotArea = {Plot_Area, Plot_Area, Plot_Area, Plot_Area};

FexNum   = {Fex_NumRoth, Fex_NumRoth, Fex_NumRoth, Fex_NumRoth};

V2Num    = {[],[],[],[]};

eqNum    = {eq_Num,eq_Num,eq_Num,eq_Num};

DDFTCode = {'DDFTDynamics', 'DDFTDynamics', 'DDFTDynamics', 'DDFTDynamics'};
        
doPlots = false;

DDFTParamsNames = {{'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','HINum','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','HINum','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','HINum','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','HINum','doPlots'}};
             
               
DDFTName={'No HI','Full HI','Only Wall', 'RP'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r','r','r','r'};

% whether to do DDFT calculations
%doDDFT={true,true,true,true}; 
doDDFT={true,true,false,true}; 
%doDDFT={false,false,false,true};

% do we load and save the DDFT data
loadDDFT={true,true,true,true};

DDFTColour = {{'r'},{'g'},{'b'},{'m'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

separateComp = true;

%viewPoint = [65;10];
viewPoint = [79;14];
%viewPoint = [0;0];

% x axis for position and velocity plots
rMin=[-10;0.5];
%rMin=[2;0];
rMax=[10;10];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.4;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[-1;6];
RMMax=[1;7];
PMMin=[-1;-1];
PMMax=[1;1];

% number of bins for histograming of stochastic data
nBins=[20;20];

% determine which movies/plots to make
% distribution movies/plots
doMovieGif     = false;          % .gif movie
doMovieAvi     = false;
doInitialFinal = false;
doMeans        = false;
doEquilibria   = false;
doSnapshotsError = true;

sendEmail = false;
