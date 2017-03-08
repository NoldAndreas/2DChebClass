%--------------------------------------------------------------------------
% compulsory physical constants
%--------------------------------------------------------------------------

geom='planar2D';

% dimension in which to do the stochastic calculations
stocDim=2;
% dimension in which to do the DDFT calculations
DDFTDim=2;

nParticlesS=[10;10];

kBT=1;          % temperature
mS=[1;1];       % particle masses

gammaS=[1;1];   % friction/diffusion
D0S=kBT./mS./gammaS;

%--------------------------------------------------------------------------
% V1 parameters
%--------------------------------------------------------------------------

V1DV1='V1_JCP_HardDiscs';

% appropriate physical parameters for potentials in V1DV1
V0S        = 0.01;
V0addS     = [3;3];
tauS       = 0.1;
sigma1AddS = 0.5;
sigma2AddS = 0.5;

y10S       = [-1;-0.75];
y20S       = [0;-0.25];
y11S       = [0.75;1];
y21S       = [0;-0.25];
y12S       = 0;
y22S       = [0.75;0.75];
% form into structure to make it easy to pass arbitrary parameters to
% potentials
potParamsNames = {'V0','V0add','tau','sigma1Add','sigma2Add',...
                  'y10','y20','y11','y21','y12','y22'};

%--------------------------------------------------------------------------
% V2 parameters
%--------------------------------------------------------------------------

V2DV2='hardSphere';

% hard sphere diameters
sigmaS = [1 1;1 1];

potParams2Names={'sigma'};

%--------------------------------------------------------------------------
% Save time setup
%--------------------------------------------------------------------------

% end time of calculation
tMax=0.25;

%--------------------------------------------------------------------------
% Stochastic setup
%--------------------------------------------------------------------------

% number of samples to take of the initial equilibrium
nSamples=5000000;  
%nSamples=5000;  
 
% intial condition for sampling
initialGuess='makeGrid';

% number of runs of stochastic dynamics to do, and average over
%nRuns=50000;
nRuns = 1;

% number of cores to use in parallel processing
%poolsize=12;
poolsize = 1;

% type of calculation, either 'rv'=Langevin or 'r'=Ermak-MCammon
stocType={'r'};

% whether to include hydrodynamic interactions
stocHI={false};
% HI interaction matrices
stocHIType={[]};

% names for stochastic calculations -- used as legend text
stocName={'r0'};

% whether to do Langevin and Brownian dynamics
doStoc={true};

% whether to load saved data for Langevin and Brownian dynamics
loadStoc={true};

% number of time steps
tSteps={10^4};

% whether to save output data (you probably should)
saveStoc={true};

% plotting colour and linestyle
stocColour = {{'g','m'}};
stocStyle = {{'--','--'}};

%--------------------------------------------------------------------------
% DDFT setup
%--------------------------------------------------------------------------

% plotting box bound
boxL = 2;

Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'L1',4,...
                   'y2Min',-inf,'y2Max',inf,'L2',4);

Plot_Area = struct('y1Min',-boxL,'y1Max',boxL,'N1',100,...
                   'y2Min',-boxL,'y2Max',boxL,'N2',100);

Fex_Num   = struct('Fex','FMTRoth',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);

eq_Num = struct('eqSolver','fsolve');
                   
paramNames = {'PhysArea','PlotArea','FexNum','V2Num','eqNum','doPlots'};

% construct computations

PhysArea = {};
PlotArea = {};
FexNum = {};
V2Num = {};
eqNum = {};
DDFTCode = {};
DDFTParamsNames = {};
DDFTType = {};
DDFTName = {};
doDDFT = {};
loadDDFT = {};

Nlist = [20;30;40;60]; % N contour plotS

%for N = 10:2:60  % error comparison
for iN = 1:length(Nlist);
    N = Nlist(iN);
    Phys_Area.N = [N,N];
    PhysArea = cat(2,PhysArea,Phys_Area);
    PlotArea = cat(2,PlotArea,Plot_Area);
    FexNum = cat(2,FexNum,Fex_Num);
    V2Num = cat(2,V2Num,{[]});
    eqNum = cat(2,eqNum,eq_Num);
    DDFTCode = cat(2,DDFTCode,'DDFTDynamics');
    DDFTParamsNames = cat(2,DDFTParamsNames,{paramNames});
    DDFTType = cat(2,DDFTType,'r');
    doDDFT = cat(2,doDDFT,true);
    loadDDFT = cat(2,loadDDFT,true);
    DDFTName = cat(2,DDFTName, num2str(N));
end

% whether to do the built-in DDFT plotting
doPlots = false;

% colour and linestyles
DDFTColour = {{'r','r'},{'g','g'},{'b','b'},{'m','m'}};
DDFTStyle = {{':',':'},{'--','--'},{'-','-'},{':',':'}};

%--------------------------------------------------------------------------
% Plotting setup
%--------------------------------------------------------------------------

plotType = 'surf';

viewPoint = [-56;7];
rMin = [-2;-1.25];
rMax = [2;1.75];

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
nBins=[50;50];

separateSpecies = true;

% determine which movies/plots to make
% distribution movies/plots
doEquilibria   = true;
