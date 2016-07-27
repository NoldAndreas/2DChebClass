function S=setDefaults(S)
% S=setDefaults(S)
%   Sets default values for parameters not in S including line styles and
%   axis limits
%
% INPUTS:
%   S  -- struct of size [1 1] containing
%              inputFile: name of input file (string)
%              tMax:      end time (scalar)
%              geom:      geometry (string)
%              anyStoc:   doing any stoc calculations (t/f)
%                 IF SO
%                 stocName: string names of stoc calculations 
%                            ({name 1, ..., name nStoc})
%                 nRuns:    number of stoc runs (scalar)  
%                 nSamples: number of initial/final samples (scalar)
%              anyDDFT:   doing any DDFT calculations (t/f)
%                 IF SO
%                 DDFTName: string names of DDFT calculations 
%                            ({name 1, ..., name DDFT})
%                 DDFTDim:  DDFT dimension (scalar)
% 
% OUTPUTS:
%   S  -- struct of size [1 1] with default parameters set for those
%          missing from S

%--------------------------------------------------------------------------
% Default plot times
%--------------------------------------------------------------------------

D.nPlots=100;
D.plotTimes=0 : S.tMax / D.nPlots : S.tMax;

%--------------------------------------------------------------------------
% Construct D0 or gamma if it's missing
%--------------------------------------------------------------------------

if(~isfield(S,'D0S'))
    D.D0S=S.kBT ./ S.mS ./ S.gammaS;
elseif(~isfield(S,'gammaS'))
    D.gammaS=S.kBT ./ S.mS ./ S.D0S;
end
    
%--------------------------------------------------------------------------
% Default stochastic options
%--------------------------------------------------------------------------

if(S.anyStoc)
    D.noise=true;
    D.MBp=true;
    D.initialGuess='makeGrid';
    D.useNewHS = false;
    D.nSamples=10*S.nRuns;
    D.thin=1;
    D.burnin=0;
    D.sampleShift=0;
    D.sampleSkip=floor((S.nSamples-D.sampleShift)/S.nRuns);
    D.poolSize=1;
    D.loadSamples=true;
    
    D.HIParamsNames={};
    D.stocUseDivergence = false;
    
    D.doStrip = false;
    
    nStoc=length(S.stocName);
    D.doStoc=cell(1,nStoc);
    D.loadStoc=cell(1,nStoc);
    D.saveStoc=cell(1,nStoc);
   
    for iStoc=1:nStoc
        D.doStoc{iStoc}=true;
        D.loadStoc{iStoc}=true;
        D.saveStoc{iStoc}=true;
    end

    if(S.stocDim==1)
        D.nBins=50;
    else
        D.nBins=[25;25];
    end       

    D.fixedBins=true;
    
    D.vCutoff=5*10^(-2);
    
end

%--------------------------------------------------------------------------
% Default DDFT options
%--------------------------------------------------------------------------

if(S.anyDDFT)
    nDDFT=length(S.DDFTName);
    D.doDDFT=cell(1,nDDFT);
    D.loadDDFT=cell(1,nDDFT);
    D.saveDDFT=cell(1,nDDFT);
    
    D.HIParamsNamesDDFT={};
    
    D.doDDFTPlots=false;
    
    D.DDFTType=cell(1,nDDFT);
   
    for iDDFT=1:nDDFT
        D.doDDFT{iDDFT}=true;
        D.loadDDFT{iDDFT}=true;
        D.saveDDFT{iDDFT}=true;
        D.DDFTType{iDDFT}='r';
    end
end

%--------------------------------------------------------------------------
% Default plotting options (all off)
%--------------------------------------------------------------------------

D.doMovieGif=false;
D.doMovieAvi=false;
D.doPdfs=false;
D.doFigs=false;
D.doMovieSwf=false;
D.doInitialFinal=false;
D.doEquilibria=false;
D.doMeans=false;
D.doDDFTPlots=false;

D.doInitialFinalP=false;
D.doMovieGifP=false;
D.doPdfsP=false;
D.doMovieSwfP=false;
D.doInitialFinalP=false;
D.doSnapshotsError = false;
D.doSnapshotsDDFT = false;

D.doCustom=false;
D.doCustomP=false;
D.custom='';

D.plotType='contour';

D.separateSpecies = false;
D.separateError = false;
D.separateComp = false;

D.legPos='off';
D.oneLeg='top';
D.perRow=4;

D.viewPoint=[-30,45];

D.symbolLabels=false;

D.dpi=300;
D.fps=5;
D.bitmap=true;
D.quiet=true;

D.plotDensity=true;
D.plotCurrent=false;

D.fixedInitial=false;
D.sampleFinal=false;
D.plotEquilibria=false;

%--------------------------------------------------------------------------
% Default line styles
%--------------------------------------------------------------------------

D.potNames=S.inputFile;
D.nParticlesS=S.nParticlesS;

if(isfield(S,'stocName'))
    D.stocName=S.stocName;
else
    D.stocName=[];
end

if(isfield(S,'DDFTName'))
    D.DDFTName=S.DDFTName;
else
    D.DDFTName=[];
end

D.geom=S.geom;

if(isfield(S,'DDFTDim'))
    D.dim=S.DDFTDim;
else
    D.dim=1;
end

D=getDefaultLines(D);

%--------------------------------------------------------------------------
% Default axis limits
%--------------------------------------------------------------------------

D=getDefaultAxes(D);


%---------

if(D.dim==1)
    D.NPlot=100;
elseif(D.dim==2)
    D.NPlot1=100;
    D.NPlot2=100;
    D.quiverSkip=10;
end

%--------------------------------------------------------------------------
% Default email
%--------------------------------------------------------------------------

D.sendEmail = false;
D.emailAddress = 'b.goddard@ed.ac.uk';

%--------------------------------------------------------------------------
% Assign values to those which don't exist in S
%--------------------------------------------------------------------------

Dfields=fields(D);
nDfields=length(Dfields);

for iField=1:nDfields
   if(~isfield(S,Dfields{iField}))
       S.(Dfields{iField})=D.(Dfields{iField});
   end
end