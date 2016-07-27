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

%inputFile = 'InfSpaceMoveN_HI_6';  % WORKS!!
%inputFile = 'HalfSpaceMoveN_HI';  % WORKS
%inputFile = 'HalfSpaceMoveAwayN_HI';
%inputFile = 'HalfSpaceMoveParallelN_HI';

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

%inputFile = 'Karolis';

%---------------------------------

%inputFile = 'PressureTest1DFree';
%inputFile = 'PressureTest1DFreeStill';
%inputFile = 'PressureTest1DFreeFade';
inputFile = 'PressureTest1DFreeFade2';
%inputFile = 'PressureTest1DGaussian';
%inputFile = 'PressureTest1DFreeHill';


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


R = stocStruct.x;
P = stocStruct.p;

R = permute(R,[2,1,3]);  %nParticles x nSamples x nTimes
P = permute(P,[2,1,3]);

nBins = 25;
optsPhys.rMin = -2;
optsPhys.rMax = 2;

for tPos = 41:10:101

    Rt = R(:,:,tPos);
    Pt = P(:,:,tPos);

    Rt = reshape(Rt,size(Rt,1),size(Rt,2)*size(Rt,3));
    Pt = reshape(Pt,size(Pt,1),size(Pt,2)*size(Pt,3));

    [PK,PV,binMids,nR,meanP,PbyBin] = getPressure1D(Rt,Pt,nBins,optsPhys);

    figure

    for iPlot = 1:nBins
        hold off
        [N,Pmid] = hist(PbyBin{iPlot},20);
        bar(Pmid,N);
        hold on
        maxN = max(N);
        plot(Pmid, maxN*exp(-(Pmid-meanP(iPlot)).^2/2) );
        title(num2str(binMids(iPlot)));
        pause
    end
    
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
