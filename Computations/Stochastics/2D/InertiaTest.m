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

gammaS=2;
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=pi*gammaS/2;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='V1_Well_osc';

% appropriate physical parameters for potentials in V1DV1
V0S        = 0.0001;
V0addS     = 2;

sigma1AddS = 5;
sigma2AddS = 5;

y10S       = 0;
y20S       = 0;

tauS       = tMax/pi;

% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'V0','V0add','tau','sigma1Add','sigma2Add',...
                  'y10','y20'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='Gaussian';

epsilonS = 2;
alphaS   = 1;

potParams2Names={'alpha','epsilon'};

%--------------------------------------------------------------------------
% HI parameters
%--------------------------------------------------------------------------

sigmaHS = 0.5;

HIParamsNames={'sigmaH'};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

Phys_Area = struct('shape','InfSpace','N',[30;30], ...
                       'y1Min',-inf,'y1Max',inf,'L1',3,...
                       'y2Min',-inf,'y2Max',inf,'L2',3);
                  
Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',-5,'y2Max',5,'N2',100);

Fex_Num = struct('Fex','Meanfield','N',[20,20],'L',1);
                   
PhysArea = {Phys_Area, Phys_Area};

PlotArea = {Plot_Area, Plot_Area};

V2Num   = {Fex_Num, Fex_Num};
Inertial = {false,true};

DDFTCode = {'DDFTDynamics', 'DDFTDynamics'};
        
doPlots = true;

DDFTParamsNames = {{'PhysArea','PlotArea','V2Num','doPlots','Inertial'}, ...
                   {'PhysArea','PlotArea','V2Num','doPlots','Inertial'}};
               

HIParamsNamesDDFT={};               
               
DDFTName={'Diffusion','Inertia'};


% type of DDFT calculations, either 'rv' to include momentum, or 'r' for
% the standard position DDFT
DDFTType={'r','r'};

% whether to do DDFT calculations
doDDFT={true,true};

% do we load and save the DDFT data
loadDDFT={true,true};
%loadDDFT={false,false};

DDFTColour = {{'r'},{'b'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

viewPoint = [-56;7];

% x axis for position and velocity plots
rMin=[-5;-5];
rMax=[5;5];
pMin=rMin;
pMax=rMax;

% y axis for position and velocity plots
RMin=0;
RMax=0.01;

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
doMovieGif     = false;          % .gif movie
doInitialFinal = true;
doMeans        = false;
doEquilibria   = false;

%sendEmail = true;
