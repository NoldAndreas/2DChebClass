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

V1DV1='V1_Well_Move_Inf_2';

% appropriate physical parameters for potentials in V1DV1
%V0S        = 0.015;
V0S        = 0.01;
%V0S        = 0.02;

nAddS = 2;

V0addS     = 3;
tauS       = 0.01;

%sigma1AddS = 1;
%sigma2AddS = 1;

sigma1AddS = 50;
sigma2AddS = 50;

y10aS       = -2;
y20aS       = -2;
y10bS       = 0;
y20bS       = 2;

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
%tMax=1;
tMax=2;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true


nSamples = 1000000;  

initialGuess='makeGrid';

sampleFinal = false;

% number of runs of stochastic dynamics to do, and average over

%nRuns = 50000;

nRuns = 10000;

%nRuns = 500;


% number of cores to use in parallel processing
poolsize=12;
%poolsize=1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r','r','r'};

% whether to include hydrodynamic interactions
stocHI={false,true,true};
stocUseDivergence = {false,true,true};


% HI interaction matrices
stocHIType={[],'RP2D_2','RP2D'};

% names for stochastic calculations -- used as legend text
stocName={'noHI','RP div 2','RP div'};

% whether to do Langevin and Brownian dynamics
%doStoc={false,false,false};
doStoc={true,false,true};

% whether to load saved data for Langevin and Brownian dynamics
%loadStoc={true,true,true};
loadStoc={true,true,true};

% number of time steps
%tSteps={10^3,10^3,10^3};
%tSteps={10^4,10^4,10^4};
tSteps={5*10^4,10^4,5*10^4};

% whether to save output data (you probably should)
saveStoc={true,true,true};

stocStyle = {{'-'},{'-'},{'-'}};
stocColour = {{'r'},{'c'},{'b'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

y1Plot=10;
y2Plot=10;

Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'N',[40,40],'L1',4,...
                    'y2Min',-inf,'y2Max',inf,'L2',4);

% Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'N',[20,20],'L1',4,...
%                     'y2Min',-inf,'y2Max',inf,'L2',4);
                
Sub_Area = struct('shape','Box','y1Min',-3,'y1Max',3,'N',[20,20],...
                      'y2Min',0.5,'y2Max',1);
                   
Plot_Area = struct('y1Min',-y1Plot,'y1Max',y1Plot,'N1',40,...
                       'y2Min',-y2Plot,'y2Max',y2Plot,'N2',40);

Fex_NumRosenfeld   = struct('Fex','FMTRosenfeld',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

Fex_NumRoth   = struct('Fex','FMTRoth',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

Fex_Num3D   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

HI_RP = struct('N',[20;20],'L',4,'HI11','noHI_2D','HI12','RP12_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D');

HIParamsNamesDDFT={'sigmaH','sigma'};                  
                  
%eq_Num    = struct('eqSolver','Newton','NewtonLambda1',0.7,'NewtonLambda2',0.7);
eq_Num = struct('eqSolver','fsolve');
                   
PhysArea = {Phys_Area, Phys_Area, Phys_Area};

SubArea  = {Sub_Area, Sub_Area, Sub_Area};

PlotArea = {Plot_Area, Plot_Area, Plot_Area};

FexNum   = {Fex_NumRosenfeld, Fex_NumRoth, Fex_NumRoth};

V2Num    = {[],[],[]};

eqNum    = {eq_Num,eq_Num,eq_Num};

HINum    = {[], ...
            [], ...
            HI_RP,...
           };

DDFTCode = {'DDFTDynamics', 'DDFTDynamics', 'DDFTDynamics'};
        
doPlots = false;

DDFTParamsNames = {{'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','HINum','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','HINum','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','HINum','doPlots'}};
             
               
DDFTName={'Rosenfeld','Roth','Roth HI'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r','r','r'};

% whether to do DDFT calculations
%doDDFT={false,true,true}; 
doDDFT={false,true,true};

% do we load and save the DDFT data
loadDDFT={true,true,true};

DDFTColour = {{'c'},{'r'},{'b'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

%separateComp = true;

%viewPoint = [65;10];
viewPoint = [-5;25];
%viewPoint = [0;0];

% x axis for position and velocity plots
rMin=[-y1Plot;-y2Plot];
%rMin=[2;0];
rMax=[y1Plot;y2Plot];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.3;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[0;-2];
RMMax=[1;2];
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
