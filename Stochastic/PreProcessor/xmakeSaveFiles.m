function fileStruct=makeSaveFiles(optsStruct,optsPhysGIF)
%fileStruct=makeSaveFiles(stocStruct,DDFTStruct,optsPhysGIF)
%   returns save file names and ensures appropriate directories exist
%   also copies input and potential files to save directory
%
% INPUTS: 
%  optsStruct -- a structure of size [1 1] containing at least
%         anyStoc       (true/false do any stoc calculations)
%         tMax          (max time)
%         tSteps        (number of time steps)
%         nRuns         (number of stochastic runs)
%         nSamples      (number of samples for intial/final distribution)
%         thin          (option for slicesampling, integer, default 1)
%         burnin        (            "           , integer, default 0)
%         stocType      (should be one of 'rv','r', string or {1,nStoc})
%         doStoc        (true/false whether to do stochastics, logical or {1,nStoc})
%         stocName      (name for stochastic calculation {1,nStoc})
%         stocHI        (true/false to include hydrodynamic interactions, logical or {1,nStoc})
%         stocHIType    (type of HI matrix to use, string or {1,nStoc})
%         stocDim       (stochastic dimension)
%         MBp           (true/false use Maxwell-Boltzmann momentum rather than zeros)
%         nBins         (number of histogramming bins (dim,1))
%
%         anyDDFT       (true/false do any stoc calculations)
%         DDFTCode      (string name of DDFT code used, string or {1,nDDFT})
%         DDFTType      (should be one of 'rv','r', string or {1,nDDFT})
%         doDDFT        (true/false whether to do DDFT calculation logical or {1,nDDFT})
%         DDFTName      (name for DDFT calculation in legends {1,nDDFT})
%         nPlots        (number of times to save plotting data)
%
%         potParamsSaveNames  (names of V1 potential parameters to save,
%                              {{'p1S',...,'pnS'},{},{}})
%         potParams2SaveNames (names of V2 potential parameters to save,
%                              {'p1S',...,'pnS'})
%         HIParamsSaveNames   (names of HI potential parameters to save,
%                              {'p1S',...,'pnS'})
%         DDFTParamsSaveNames (names of DDFT parameters to save, either
%                              {{'p1',...,'pn'}} or 
%                              {{'p1',...,'pn'}, ..., {...}} )
%         DDFTParams          (struct 1xnDDFT, contents as in
%                               DDFTParamsSaveNames)
%         V1DV1               (V1 potentials, {V1G, V1I, V1F})
%         V2DV2               (V2 potentials, string)
%         inputFile           (input file name, string)
%
%  optsPhysGIF -- a structure of size [1 3] containing at least:
%         tMax          (max time)
%         geom          (geometry; string')
%         dim           (dimension)
%         D0S           (diffusion constant, (nSpecies,1))
%         gammaS        (friction constant, (nSpecies,1))
%         kBT           (temperature)
%         nParticlesS   (number of particle in each species, (nSpecies,1))
%         mS            (particle mass, (nSpecies,1))
%         V1DV1         (strings identifying potential and gradient functions
%                        {'General','Initial','Final'})
%         potNames      (save names for potentials)
%         <V1 parameters>     (as in potParamsSaveNames (nSpecies,1))
%         <V2 parameters>     (as in potParams2SaveNames (nSpecies,nSpecies))
%         <HI parameters>     (as in HIParamsSaveNames (nSpecies,1))
%
% OUTPUTS:
%  fileStruct -- a structure of size [1 1] containing
%    saveFileI     --  a string location for the initial sampling data file
%    saveFileF     --  a string location for the final sampling data file
%    saveFilePIF   --  a string location for the momentum sampling data file
%    saveFileIEq   --  a string location for the initial equilibria data file
%    saveFileFEq   --  a string location for the final equilibria data file
%    saveFileStoc  --  a nStoc long cell of string locations for the 
%                      stochastic dynamics data files
%    saveFileDDFT  --  a nDDFT long cell of string locations for the DDFT 
%                      dynamics data files
%    plotFile      --  a length 3 cell of string locations for the  
%                       (1) general, (2) intial and (3) final plots
%    plotFileP     --  a length 3 cell of string locations for the  
%                       (1) general, (2) intial and (3) final particle plots
%    movieFile     --  a string location for the movie file.  No extension!
%    movieFileP    --  a string location for the particle movie file.  No extension!
%    pdfDir        --  a string directory name for saving pdfs to make swf
%    pdfDirP       --  a string directory name for saving particle pdfs to 
%                       make swf
%    saveFile      --  placeholder, []

%----------------------------------------------------------------------
%  Get relevant physical information
%----------------------------------------------------------------------

% description of the 1-, 2-body and HI 

potDescrGIF=cell(1,3);

% possibly different for General, Initial and Final potentials (only in the
% stochastic case - DDFT calculations use one time-dependent potential

anyStoc = optsStruct.anyStoc;

for iDescr=1:3

    % save name of potential
    potDescrGIF{iDescr}=optsPhysGIF(iDescr).potNames;
    
    % V1 parameters to include in save name
    potParamsNames=optsStruct.potParamsSaveNames{iDescr};
    
    % construct list of parameter values to include in save name
    nParams=length(potParamsNames);
    for iParam=1:nParams
        % name of this parameter
        potName=potParamsNames{iParam};
        % and its value
        if(anyStoc)
            values=optsPhysGIF(iDescr).(potName);
        else
            values=optsPhysGIF(iDescr).V1.(potName);
        end
        % for multiple species there's more than one value for each
        % parameter
        nValues=length(values);
        for iValue=1:nValues
            potDescrGIF{iDescr}=cat(2,potDescrGIF{iDescr},['-' num2str(values(iValue))]);
        end

    end    
    
    % V2 parameters to include in save name
    potParams2Names=optsStruct.potParams2SaveNames;
    
    % construct list of parameter values to include in save name
    nParams2=length(potParams2Names);
    for iParam=1:nParams2
        potName=potParams2Names{iParam};
        if(anyStoc)
            values=optsPhysGIF(iDescr).(potName);
        else
            values=optsPhysGIF(iDescr).V2.(potName);
        end
        values=values(:);
        nValues=length(values);
        for iValue=1:nValues
            potDescrGIF{iDescr}=cat(2,potDescrGIF{iDescr},['-' num2str(values(iValue))]);
        end

    end
    
    % HI parameters to include in save name
    HIParamsNames=optsStruct.HIParamsSaveNames;

    % construct list of parameter values to include in save name
    nParams=length(HIParamsNames);
    for iParam=1:nParams
        HIName=HIParamsNames{iParam};
        if(anyStoc)
            values=optsPhysGIF(iDescr).(HIName);
        else
            values=optsPhysGIF(iDescr).HI.(HIName);
        end
        nValues=length(values);
        for iValue=1:nValues
            potDescrGIF{iDescr}=cat(2,potDescrGIF{iDescr},['-' num2str(values(iValue))]);
        end
    end    
    
end

% if all three descriptions are equal (e.g. in a DDFT computation) then
% only save the name once
if( isequal(potDescrGIF{1},potDescrGIF{2}) && isequal(potDescrGIF{1},potDescrGIF{3}) )

    potDescr=[ potDescrGIF{1} ];
    
else
    
    % otherwise save the full description
    potDescr=[ potDescrGIF{1} potDescrGIF{2} potDescrGIF{3}];
    
end


% get compulsory physical parameters for naming save directory
nParticlesS=optsPhysGIF(1).nParticlesS;
mS=optsPhysGIF(1).mS;
gammaS=optsPhysGIF(1).gammaS;
D0S=optsPhysGIF(1).D0S;

kBT=optsPhysGIF(1).kBT;
geom=optsPhysGIF(1).geom;
tMax=optsPhysGIF(1).tMax;
stocDim=optsPhysGIF(1).dim;

nParticlesText=[];
mText=[];
gammaText=[];
D0Text=[];

% construct description for all species
nSpecies=length(nParticlesS);
for iSpecies=1:nSpecies
    nParticlesText=cat(2,nParticlesText,['-' num2str(nParticlesS(iSpecies))] );
    mText=cat(2,mText,['-' num2str(mS(iSpecies))] );
    gammaText=cat(2,gammaText,['-' num2str(gammaS(iSpecies))] );
    D0Text=cat(2,D0Text,['-' num2str(D0S(iSpecies))] );
end

% standard  physical description
physDescr=['-' geom nParticlesText mText gammaText D0Text ...
           '-' num2str(tMax) '-' num2str(kBT)  ...
           '-' num2str(stocDim)];

%----------------------------------------------------------------------
% Set up save file directories
%----------------------------------------------------------------------

% stochastic information
if(optsStruct.anyStoc)
    stocStruct=struct('tMax',optsStruct.tMax,'tSteps',optsStruct.tSteps,'nRuns',optsStruct.nRuns, ...
                'nSamples',optsStruct.nSamples,'sampleShift',optsStruct.sampleShift, ...
                'sampleSkip',optsStruct.sampleSkip,'thin',optsStruct.thin,'burnin',optsStruct.burnin, ...
                'stocType',optsStruct.stocType,'doStoc',optsStruct.doStoc,...
                'stocName',optsStruct.stocName,'stocHI',optsStruct.stocHI, 'stocHIType', optsStruct.stocHIType, ...
                'stocDim',optsStruct.stocDim,'nPlots',optsStruct.nPlots,'MBp',optsStruct.MBp,'nBins',optsStruct.nBins);
else
    stocStruct=[];
end
            
% DDFT information
if(optsStruct.anyDDFT)            
    DDFTStruct=struct('tMax',optsStruct.tMax,'DDFTParamsSaveNames',optsStruct.DDFTParamsSaveNames, ...
                'DDFTType',optsStruct.DDFTType,'DDFTDir',optsStruct.DDFTDir,'DDFTCode',optsStruct.DDFTCode, ...
                 'doDDFT',optsStruct.doDDFT,'DDFTName',optsStruct.DDFTName,'nPlots',optsStruct.nPlots);            
else
    DDFTStruct=[];
end


% number of stochastic calculations
nStoc=size(stocStruct,2);
% number of DDFT calculations
nDDFT=size(DDFTStruct,2);

% create main save directory
saveDir=['Data' filesep potDescr physDescr filesep]; 
if(~exist(saveDir,'dir'))
    mkdir(saveDir);
end

if(nStoc>0)
    % create sampling directory
    sampleDir=['Data' filesep 'Equilibrium' filesep];
    if(~exist(sampleDir,'dir'))
        mkdir(sampleDir);
    end

    % and a place to save the corresponding means
    sampleMeanDir=['Data' filesep 'Equilibrium' filesep 'Means' filesep];
    if(~exist(sampleMeanDir,'dir'))
        mkdir(sampleMeanDir);
    end

    % and one for momentum sampling 
    sampleMomDir=['Data' filesep 'Equilibrium' filesep 'Momentum' filesep];
    if(~exist(sampleMomDir,'dir'))
        mkdir(sampleMomDir);
    end

    % create stochastic data directory
    stocDir=[saveDir 'Stochastic' filesep];
    if(~exist(stocDir,'dir'))
        mkdir(stocDir);
    end
    
    % create stochastic mean directory
    stocMeanDir=[stocDir 'Means' filesep];
    if(~exist(stocMeanDir,'dir'))
        mkdir(stocMeanDir);
    end
    
end
    
if(nDDFT>0)    
    % create DDFT data directory
    DDFTDir=[saveDir 'DDFT' filesep];
    if(~exist(DDFTDir,'dir'))
        mkdir(DDFTDir);
    end
end

% create output (plotting) directory
outputDir=[saveDir 'Output' filesep];
if(~exist(outputDir,'dir'))
    mkdir(outputDir);
end 
       
%----------------------------------------------------------------------
% Set up stochastic equilibrium save files
%----------------------------------------------------------------------

if(nStoc>0)
    % whether to do Maxwell-Boltzmann or set momenta to zero
    MBp=stocStruct.MBp;
    
    % stochastic sampling options
    nSamples=stocStruct.nSamples;
    thin=stocStruct.thin;
    burnin=stocStruct.burnin;
    nBins=stocStruct.nBins;
    
    % number of bins text
    nBinsText=[];
    
    for iDim=1:length(nBins)
        nBinsText=cat(2,nBinsText,['-' num2str(nBins(iDim))] );
    end
    
    % data relevant to sampling
    sampleDescr=[num2str(nSamples) '-' num2str(thin) '-' num2str(burnin) ...
                 nParticlesText mText '-' num2str(kBT) ];

    % for initial and final sampling data, note it doesn't actually matter if
    % we sampled as an initial or a final case
    saveFileI=[sampleDir potDescrGIF{2} '-0-' sampleDescr '.mat'];
    saveFileF=[sampleDir potDescrGIF{3} '-' num2str(tMax) '-' sampleDescr '.mat'];

    % equilibrium data, include number of bins as this affects equilibrium
    % plots
    saveFileIEq=[sampleMeanDir potDescrGIF{2} '-0-Eq-' sampleDescr '-' nBinsText '.mat'];
    saveFileFEq=[sampleMeanDir potDescrGIF{3} '-' num2str(tMax) '-Eq-' sampleDescr '-' nBinsText '.mat'];
    
    % sampling for momentum
    % data relevant to momentum sampling (note do not need sigma)
    sampleDescrMom=['-' num2str(nSamples) '-' num2str(thin) '-' num2str(burnin) ...
                 nParticlesText mText '-' num2str(kBT) '-' num2str(stocDim) ];
    saveFilePIF=[sampleMomDir 'momentum-' num2str(MBp) '-' sampleDescrMom '.mat'];
end
    
    
%----------------------------------------------------------------------
% Set up stochastic descriptions and save files
%----------------------------------------------------------------------

if(nStoc>0)
    % preallocate for naming 
    saveFileStoc=cell(1,nStoc);
    saveFileStocMeans=cell(1,nStoc);

    % for naming output file
    stocDescrFull=[ ];

    for iStoc=1:nStoc

        % stochastic dynamics options
        stocType=stocStruct(iStoc).stocType;
        stocHI=stocStruct(iStoc).stocHI;
        tSteps=stocStruct(iStoc).tSteps;
        nRuns=stocStruct(iStoc).nRuns;
        nPlots=stocStruct(iStoc).nPlots;
        % name for HI
        if(stocHI)
            stocHIType=stocStruct(iStoc).stocHIType;
        else
            stocHIType='Off';
        end

        % save file stochastic information
        stocDescr=[stocType '-' num2str(MBp) '-' stocHIType '-' num2str(tSteps)...
                           '-' num2str(nRuns) '-' num2str(nPlots)];

        % file in which to save stochastic data
        saveFileStoc{iStoc}=[stocDir stocDescr '.mat']; 
        saveFileStocMeans{iStoc}=[stocMeanDir stocDescr '.mat']; 

        % if we're going to do the calculation, add it to the string for the
        % output file
        if(stocStruct(iStoc).doStoc)
           stocDescrFull=cat(2,stocDescrFull, '-'); 
           stocDescrFull=cat(2,stocDescrFull, stocDescr);
        end

    end
end

%----------------------------------------------------------------------
% Set up DDFT descriptions and save files
%----------------------------------------------------------------------

if(nDDFT>0)
    % preallocate for individual files
    saveFileDDFT=cell(1,nDDFT);

    % for naming output file
    DDFTDescrFull=[ ];

    for iDDFT=1:nDDFT

        % DDFT options
        DDFTCode=DDFTStruct(iDDFT).DDFTCode;
        DDFTType=DDFTStruct(iDDFT).DDFTType;
        nPlots=DDFTStruct(iDDFT).nPlots; 
        
        DDFTParamsSaveNames=DDFTStruct(iDDFT).DDFTParamsSaveNames;
        nParams=length(DDFTParamsSaveNames);
        DDFTParams=optsStruct.DDFTParams(iDDFT);

        % DDFT parameter description
        parDescr=[];
        for iParam=1:nParams
            
            % names and values of parameters
            parName=DDFTParamsSaveNames{iParam};
            values=DDFTParams.(parName);
            
            if(ischar(values))
                % if it's a string (e.g. HI name) then save the string
                parDescr=cat(2,parDescr,['-' values]);   
            else
                % othereise save the values
                nValues=length(values);
                for iValue=1:nValues   
                    parDescr=cat(2,parDescr,['-' num2str(values(iValue))]);
                end
            end
        end    
        % save file DDFT information
        DDFTDescr=[DDFTCode '-' DDFTType  ...
                          parDescr '-' num2str(nPlots)];

        % file in which to save DDFT data
        saveFileDDFT{iDDFT}=[DDFTDir DDFTDescr '.mat']; 

        % if we're going to do the calculation, add it to the string for the
        % output file
        if(DDFTStruct(iDDFT).doDDFT)
           DDFTDescrFull=cat(2,DDFTDescrFull, DDFTDescr);
        end

    end
end

%----------------------------------------------------------------------
% Set up plotting save files
%----------------------------------------------------------------------

% note no extension, add this in when determining save type
movieFile=[outputDir 'dynamics'];
movieFileP=[outputDir 'particleDynamics'];
plotFile={[outputDir 'means.pdf'], ...
        [outputDir 'initial.pdf'], ...
        [outputDir 'final.pdf']};
plotFileP={[outputDir 'particleMeans.gif'], ...
    [outputDir 'particleIF.gif'], ...
    [outputDir 'particleIF.gif']}; 

% create pdf data directories
pdfDir=[saveDir 'Output' filesep 'pdfs' filesep];
pdfDirP=[saveDir 'Output' filesep 'pdfsP' filesep];

%----------------------------------------------------------------------
% Form output structure
%----------------------------------------------------------------------

% stochastic files
if(nStoc>0)
    fileStruct.saveFileI=saveFileI;
    fileStruct.saveFileF=saveFileF;
    fileStruct.saveFilePIF=saveFilePIF;
    fileStruct.saveFileIEq=saveFileIEq;
    fileStruct.saveFileFEq=saveFileFEq;
    fileStruct.saveFileStoc=saveFileStoc;
    fileStruct.saveFileStocMeans=saveFileStocMeans;
end

% DDFT files
if(nDDFT>0)
    fileStruct.saveFileDDFT=saveFileDDFT;
end

% general (e.g. plotting) files
fileStruct.saveDir=saveDir;
fileStruct.plotFile=plotFile;
fileStruct.movieFile=movieFile;
fileStruct.plotFileP=plotFileP;
fileStruct.movieFileP=movieFileP;
fileStruct.pdfDir=pdfDir;
fileStruct.pdfDirP=pdfDirP;
fileStruct.saveFile=[];

%----------------------------------------------------------------------
% Copy files to data directory
%----------------------------------------------------------------------

copyfile(['Inputs/' optsStruct.inputFile '.m'],[fileStruct.saveDir '/input.m'])
copyfile(['Potentials/' optsStruct.V1DV1{1} '.m'],[fileStruct.saveDir '/V1DV1G.m'])
copyfile(['Potentials/' optsStruct.V1DV1{2} '.m'],[fileStruct.saveDir '/V1DV1I.m'])
copyfile(['Potentials/' optsStruct.V1DV1{3} '.m'],[fileStruct.saveDir '/V1DV1F.m'])
copyfile(['Potentials/' optsStruct.V2DV2 '.m'],[fileStruct.saveDir '/V2DV2.m'])

