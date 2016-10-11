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

V1DV1='V1_Well_Move_Inf';

% appropriate physical parameters for potentials in V1DV1
V0S        = 0.01;
%V0addS     = 1.5;
V0addS     = 2;
tauS       = 0.01;

sigma1AddS = 1;
sigma2AddS = 1;

y10aS       = 2;
y20aS       = 2;
y10bS       = 0;
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

sigmaHS = 1;

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

%burnin = 20000;
%nSamples = 1000000;  
 
nSamples = 1000000;  

initialGuess='makeGridPos';

sampleFinal = false;

% number of runs of stochastic dynamics to do, and average over

%nRuns = 50000;

nRuns = 5000;

% number of cores to use in parallel processing
poolsize=16;
%poolsize=1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r','r','r'};

% whether to include hydrodynamic interactions
stocHI={false,true,true};
% HI interaction matrices
stocHIType={[],'RP','OseenPlusWall2D'};

% names for stochastic calculations -- used as legend text
stocName={'noHI','RP','OseenWall'};

% whether to do Langevin and Brownian dynamics
%doStoc={true,true,false};
doStoc={false,false,false};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true,true,true};

% number of time steps
%tSteps={10^4,10^3,10^3};
tSteps={10^4,10^4,10^3};

% whether to save output data (you probably should)
saveStoc={true,true,true};

stocStyle = {{'-'},{'-'},{'-'}};
stocColour = {{'g'},{'m'},{'b'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'N',[30,30],'L1',4,...
                   'y2Min',-inf,'y2Max',inf,'L2',4);

Sub_Area = struct('shape','Box','y1Min',-3,'y1Max',3,'N',[20,20],...
                      'y2Min',0.5,'y2Max',1);
                   
Plot_Area = struct('y1Min',-10,'y1Max',10,'N1',50,...
                       'y2Min',-10,'y2Max',10,'N2',50);

Fex_NumRosenfeld   = struct('Fex','FMTRosenfeld',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

Fex_NumRoth   = struct('Fex','FMTRoth',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

Fex_Num3D   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

HI_RP = struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','RP12_2D', ...
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
doDDFT={false,true,true}; 
%doDDFT={false,true,false};

% do we load and save the DDFT data
loadDDFT={true,true,true};

DDFTColour = {{'r'},{'b'},{'r'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

viewPoint = [-56;7];

% x axis for position and velocity plots
rMin=[-3;0];
%rMin=[2;0];
rMax=[3;3];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.4;

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
doMovieAvi     = true;
doInitialFinal = false;
doMeans        = false;
doEquilibria   = false;

%sendEmail = true;
