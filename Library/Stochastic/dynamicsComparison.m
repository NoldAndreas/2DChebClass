clear all

%--------------------------------------------------------------------------
% Set up necessary paths
%--------------------------------------------------------------------------

AddPaths();

%--------------------------------------------------------------------------
% Choose input file
%--------------------------------------------------------------------------

% CONFINED HI PAPER
inputFile = 'ConfinedHIPaperInfSpace';  
%inputFile = 'ConfinedHIPaperHalfSpaceTowards';  
%inputFile = 'ConfinedHIPaperHalfSpaceAway'; 

% COMPUTATIONAL PAPER
%inputFile = 'JCP_639_2017_Fig10_14_Gaussians';
%inputFile = 'JCP_639_2017_Fig11_15_HardDiscs';

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
