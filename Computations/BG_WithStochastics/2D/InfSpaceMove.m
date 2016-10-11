%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations (=1), but we need to know
% it's one in certain places
DDFTDim=2;

%nParticlesS=20;
nParticlesS=10;

kBT=1;          % temperature
mS=1;

gammaS=1;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='V1_Well_Move_Inf';

% appropriate physical parameters for potentials in V1DV1
V0S        = 0.01;
V0addS     = 3;
tauS       = 0.1;
%sigma1AddS = 0.5;
%sigma2AddS = 0.5;

sigma1AddS = 1;
sigma2AddS = 1;

% y10aS       = -2;
% y20aS       = 2;
% y10bS       = -2;
% y20bS       = 3;

% y10aS       = -2;
% y20aS       = 2;
% y10bS       = 2;
% y20bS       = 2;

y10aS       = 2;
y20aS       = 2;
y10bS       = -2;
y20bS       = 2;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'V0','V0add','tau','sigma1Add','sigma2Add',...
                  'y10a','y20a','y10b','y20b'};

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
tMax=0.5;


%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial and final equilibrium
% distributions goverened by the second and third arguments of V1DV1 above
% only relevant if fixedInitial=false or sampleFinal=true

%burnin = 20000;
nSamples = 1000000;  
 

%nSamples = 50000;

initialGuess='makeGridPos';

sampleFinal = false;

% number of runs of stochastic dynamics to do, and average over

nRuns=50000;

%nRuns = 1000;

% number of cores to use in parallel processing
poolsize=16;
%poolsize=1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r','r','r'};

% whether to include hydrodynamic interactions
stocHI={false,false,true};
% HI interaction matrices
stocHIType={[],'RP','OseenPlusWall2D'};

% names for stochastic calculations -- used as legend text
stocName={'noHI','RP','OseenWall'};

% whether to do Langevin and Brownian dynamics
doStoc={false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true};

% number of time steps
%tSteps={10^4,10^3,10^3};
tSteps={10^4,10^3,10^3};

% whether to save output data (you probably should)
saveStoc={true,true,true};

stocColour = {{'g'},{'g'},{'b'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'N',[30,30],'L1',4,...
                   'y2Min',-inf,'y2Max',inf,'L2',4);

Sub_Area = struct('shape','Box','y1Min',-3,'y1Max',3,'N',[20,20],...
                      'y2Min',0.5,'y2Max',1);
                   
Plot_Area = struct('y1Min',-3,'y1Max',3,'N1',100,...
                       'y2Min',0.5,'y2Max',20,'N2',100);

Fex_Num1   = struct('Fex','FMTRosenfeld',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

Fex_Num2   = struct('Fex','FMTRoth',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

Fex_Num3   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

%eq_Num    = struct('eqSolver','Newton','NewtonLambda1',0.7,'NewtonLambda2',0.7);
eq_Num = struct('eqSolver','fsolve');
                   
PhysArea = {Phys_Area, Phys_Area, Phys_Area};

SubArea  = {Sub_Area, Sub_Area, Sub_Area};

PlotArea = {Plot_Area, Plot_Area, Plot_Area};

FexNum   = {Fex_Num1, Fex_Num2, Fex_Num3};

V2Num    = {[],[],[]};

eqNum    = {eq_Num,eq_Num,eq_Num};

HINum    = {[], [],[]};

DDFTCode = {'DDFTDynamics', 'DDFTDynamics', 'DDFTDynamics'};
        
doPlots = false;

DDFTParamsNames = {{'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','doPlots'}, ...
                   {'PhysArea','SubArea','PlotArea','FexNum','V2Num','eqNum','doPlots'}};

HIParamsNamesDDFT={};               
               
DDFTName={'Rosenfeld','Roth','Rosenfeld3D'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r','r','r'};

% whether to do DDFT calculations
doDDFT={true,true,false};   % 3D case doesn't seem to work?
%doDDFT={false,true,false};

% do we load and save the DDFT data
loadDDFT={true,true,false};

DDFTColour = {{'r'},{'b'},{'g'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

viewPoint = [-56;7];

% x axis for position and velocity plots
rMin=[-3;0];
%rMax=[3;20];
rMax=[3;5];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.6;

PMin=[-1;-1];
PMax=[1;1];

% y axis for mean position and velocity plots
RMMin=[-3;5];
RMMax=[3;7];
PMMin=[-1;-1];
PMMax=[1;1];

% number of bins for histograming of stochastic data
nBins=[40;40];

% determine which movies/plots to make
% distribution movies/plots
doMovieGif     = false;          % .gif movie
doMovieAvi     = true;
doInitialFinal = false;
doMeans        = false;
doEquilibria   = false;

%sendEmail = true;
