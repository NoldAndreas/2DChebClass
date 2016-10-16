clear all

%--------------------------------------------------------------------------
% Set up necessary paths
%--------------------------------------------------------------------------

AddPaths();

%--------------------------------------------------------------------------
% Choose input file
%--------------------------------------------------------------------------

%inputFile = 'sedimentation9New';

%inputFile='APS12HS50new';
%inputFile = 'Eric300';

%inputFile='APS12HS50';
%inputFile='APS12G50';

%inputFile = 'ThreeSpeciesTest';
%inputFile = 'ThreeSpeciesTestEq';

%inputFile = 'HIWallTestTowards';  % CHECK
%inputFile = 'HIWallTestAway';
%inputFile = 'HIWallTestTowardsUnbounded';  % Shows large error when using
                                            % Rosenfeld3D
%inputFile = 'HIWallTestTowardsFMTTest';  % CHECK

%-----------------------------------

% Validation of 2D dynamics

%inputFile = 'HalfSpaceMove'; % WORKS!
%inputFile = 'InfSpaceMove'; % WORKS!

%inputFile = 'HalfSpaceMoveN'; 
%inputFile = 'InfSpaceMoveN'; % WORKS! - needs more samples?
%inputFile = 'InfSpaceMoveN_HI';
%inputFile = 'InfSpaceMoveN_HI_2';

%inputFile = 'InfSpaceMoveN_HI_3';
%inputFile = 'InfSpaceMoveN_HI_4';
%inputFile = 'InfSpaceMoveN_HI_5';

%inputFile = 'HalfSpaceTestHIDiv';

%inputFile = 'InfSpaceMoveN_HI_6_TestTiming';
%inputFile = 'HalfSpaceMoveN_HI_Short';


% test Newton vs fsolve
%inputFile = 'InfSpaceTestEq';
%inputFile = 'HalfSpaceTestEq';

%-----------------------------------

%inputFile = 'HIWallTestTowardsMove';
%inputFile = 'HIWallTestAwayMove';  
%inputFile = 'HIWallTestTowardsMoveNew';


%inputFile = 'InertiaTest';
%inputFile = 'InertiaTestInfDisc';

%inputFile = 'HITest';
%inputFile = 'HITest_2Species';

%inputFile = 'GaussianTest';
%inputFile = 'GaussianTest_2Species';

%inputFile = 'GaussianBoxTest';
%inputFile = 'GaussianBoxTest2';
%inputFile = 'GaussianBoxTest3';

%inputFile = 'BoxTest3'; % for code paper
%inputFile = 'BoxTest3N';

%inputFile = 'FMTTest_2Species'; % higher density, less accurate
%inputFile = 'FMTTest_2Species2';
%inputFile = 'FMTTest_2Species2N';

%inputFile = 'BoxTestFlow';
%inputFile = 'BoxTestFlowIdeal';  % Want to check with large scale dynamics

%inputFile = 'FMTTest_Unbounded'; % Want to check with large scale dynamics

%inputFile = 'FlowCheck';

%---------------------------------

%inputFile = 'PressureTest1DFree';
%inputFile = 'PressureTest1DFreeStill';
%inputFile = 'PressureTest1DFreeFade';
%inputFile = 'PressureTest1DFreeFade2';  % Good example for failure of LE
%inputFile = 'PressureTest1DGaussian';
%inputFile = 'PressureTest1DFreeHill';
%inputFile = 'PressureTest1DGaussianHill';

%---------------------------------

%inputFile = 'Karolis';
%inputFile = 'Box3DSampling';

%---------------------------------

% CONFINED HI PAPER

%inputFile = 'InfSpaceMoveN_HI_6';  % For confined HI paper
%inputFile = 'HalfSpaceMoveN_HI';  % For confined HI paper
%inputFile = 'HalfSpaceMoveAwayN_HI'; % For confined HI paper
%inputFile = 'HalfSpaceMoveParallelN_HI';

%inputFile = 'ConfinedHIPaperInfSpace';  % For confined HI paper
%inputFile = 'ConfinedHIPaperHalfSpaceTowards';  % For confined HI paper
inputFile = 'ConfinedHIPaperHalfSpaceAway'; % For confined HI paper

%--------------------------------------------------------------------------
% Get parameters from input file
%--------------------------------------------------------------------------

eval(inputFile);

%--------------------------------------------------------------------------
% Combine all variables in workspace into a structure for easy access
%--------------------------------------------------------------------------

vars = evalin('base','whos');
for iVar = 1:length(vars)
    optsStruct.(vars(iVar).name) = evalin('base',vars(iVar).name);
end

%--------------------------------------------------------------------------
% Run the preprocessor
%--------------------------------------------------------------------------

optsStruct=preProcess(optsStruct);

if (~optsStruct.anyStoc && ~optsStruct.anyDDFT)
    error('dynamicsComparison:noCalculations', ...
            'No calculations to be run!')
end

%--------------------------------------------------------------------------
% Collate physical information
%--------------------------------------------------------------------------

[optsPhys,optsPhysDDFT] = getOptsPhys(optsStruct);

%--------------------------------------------------------------------------             
% Collate options for stochastic and DDFT calculations
%--------------------------------------------------------------------------

optsStoc = getOptsStoc(optsStruct);

optsNumDDFT = getOptsDDFT(optsStruct);

%--------------------------------------------------------------------------
% Collate options for plotting
%--------------------------------------------------------------------------
   
[optsPlot,optsPlotParticles] = getOptsPlot(optsStruct);

%--------------------------------------------------------------------------
% Stochastic calculation
%--------------------------------------------------------------------------

if(~isempty(optsStoc))  
    stocStruct = doStochastics(optsPhys,optsStoc);
else
    stocStruct = struct([]);
end

%--------------------------------------------------------------------------
% DDFT calculation
%--------------------------------------------------------------------------

if(~isempty(optsNumDDFT))
    DDFTStruct = doDDFTs(optsNumDDFT,optsPhysDDFT,optsPlot);
else
    DDFTStruct = struct([]);
end

%--------------------------------------------------------------------------
% File list storage
%--------------------------------------------------------------------------

saveFileList(inputFile,optsPhys.potNames,stocStruct,DDFTStruct);

%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------

% in case doPlots is a variable from anything
if(exist('doPlots','var'))    
    clear('doPlots');
end

if(optsStruct.anyPlots || optsStruct.anyPlotsP)
   plotFiles = doPlots(stocStruct,DDFTStruct,optsStoc,optsNumDDFT,optsPlot,optsPlotParticles,optsPhys);
end

%--------------------------------------------------------------------------
% Pressure tensor testing
%--------------------------------------------------------------------------

%doLE = true;
doLE = false;

if(doLE)

    testLE(stocStruct,optsPhys);
    
end
    
%--------------------------------------------------------------------------
% Send Email
%--------------------------------------------------------------------------

if(optsStruct.sendEmail)
    setupEmail();
    if(~isempty(plotFiles))
        subject = ['Code complete: ' inputFile];
        body    = 'Files attached.';
        sendmail(optsStruct.emailAddress, subject, body, plotFiles);
    else
        subject = ['Code complete: ' inputFile];
        body    = 'No plot files produced.';
        sendmail(optsStruct.emailAddress, subject, body);
    end
end
