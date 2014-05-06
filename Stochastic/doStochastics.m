function stocStruct=doStochastics(optsPhys,optsStocFull)
% [stocStruct,xInitial,xFinal,pIF]=doStochastics(optsPhysGIF,optsStocVR,optsPlotGIF) 
%   determines which sampling and stochastic dynamics calculations to
%   perform and does so.
%
% INPUTS: 
%  optsPhysGIF -- a structure of size [1 3] containing at least:
%         geom          (geometry; 'spherical'/'planar')
%         dim           (dimension; 1 or 3)
%         D0            (diffusion constant)
%         kBT           (temperature)
%         nParticles    (number of particles)
%         sigma         (particle diameter)
%         m             (particle mass)
%         V1DV1         (strings identifying potential and gradient functions
%                        {'General','Initial','Final'})
%         [pot params]  (parameters for potential where
%                           optsPhysGIF(1)=dynamic
%                           optsPhysGIF(2)=initial
%                           optsPhysGIF(3)=final)
%
% optsStocFull -- a structure of size [1 nStoc] containing at least
%         tMax          (max time)
%         tSteps        (number of time steps)
%         nRuns         (number of stochastic runs)
%         nSamples      (number of samples for intial/final distribution)
%         sampleShift   (runIC(i) = sample(sampleShift+sampleSkip*i))
%         sampleSkip    (                    "                      )
%         poolsize      (number of parallel processors to use)
%         thin          (option for slicesampling, integer, default 1)
%         burnin        (            "           , integer, default 0)
%         saveFileI     (string, initial distribution save file)
%         saveFileF     (string, final distribution save file)
%         saveFile      (string, temporary save file)
%         saveFileData  (string, dynamics save file)
%         saveData      (true/false, whether to save dynamics data)
%         loadData      (true/false, whether to load dynamics data)
%         loadSamples   (true/false, whether to load sampling data)
%         initialGuess  (string, initial guess function)
%         type          (should be one of 'rv','r')
%         HI            (true/false to include hydrodynamic interactions)
%         doStoc        (true/false whether to do stochastics)
%         fixedInitial  (true/false, fixed IC or sample via potential)
%         sampleFinal   (true/false, sample final distribution or not)
%         plotTimes     (vector of save/plotting times)
%       
% OUTPUTS:
%  stocStruct -- a structure of size [1 nStoc] containing
%           x           -- a matrix of positions of size 
%                             (nRuns, dim x nParticles, length(plotTimes))
%           p           -- a matrix of momenta of size 
%                             (nRuns, dim x nParticles, length(plotTimes))
%  xInitial    -- a matrix of sampled initial positions of size 
%                  (nSamples, dim x nParticles)
%  xFinal      -- a matrix of sampled final positions of size 
%                  (nSamples, dim x nParticles)
%  pIF         -- a matrix of sampled momenta of size 
%                  (nSamples, dim x nParticles)


% number of stochastic calculations to do
nStoc=length(optsStocFull);

% determine whether we're doing any momentum calculations
doP=false;
for iStoc=1:nStoc
    if(strcmp(optsStocFull(iStoc).type,'rv'))
        doP=true;
    end
end

% set options to be those for the first stochastic calculation.  Note that 
% the values we use (nRuns, nSamples) should be the same for all
optsStoc=optsStocFull(1);

sampleFinal = optsStoc.sampleFinal;
loadSamples = optsStoc.loadSamples;
MBp         = optsStoc.MBp;


% optsStoc = rmfield(optsStoc,'loadSamples');
% optsStoc = rmfield(optsStoc,'fixedInitial');
% optsStoc = rmfield(optsStoc,'sampleFinal');
% optsStoc = rmfield(optsStoc,'poolsize');
% optsStoc = rmfield(optsStoc,'doStoc');
% optsStoc = rmfield(optsStoc,'loadStoc');
% optsStoc = rmfield(optsStoc,'saveStoc');
% optsStoc = rmfield(optsStoc,'HI');
% optsStoc = rmfield(optsStoc,'HIType');
% optsStoc = rmfield(optsStoc,'stocName');
% optsStoc = rmfield(optsStoc,'noise');
% optsStoc = rmfield(optsStoc,'nRuns');
% optsStoc = rmfield(optsStoc,'nRuns');

optsStocI.initialGuess = optsStoc.initialGuess;
optsStocI.nSamples     = optsStoc.nSamples;
optsStocI.thin         = optsStoc.thin;
optsStocI.burnin       = optsStoc.burnin;

%--------------------------------------------------------------------------
% Initial distribution
%--------------------------------------------------------------------------

fprintf(1,'Starting initial sampling ... ');

optsPhys.t=0;
opts.optsPhys = optsPhys;
opts.optsStoc = optsStocI;

ICDir = [optsPhys.potNames filesep 'Stochastic' filesep 'Initial'];

xInitial = DataStorage(ICDir,@samplepdf,opts,[],~loadSamples);

fprintf(1,'Finished\n');

%--------------------------------------------------------------------------
% Final distribution
%--------------------------------------------------------------------------

if(sampleFinal)
    fprintf(1,'Starting final sampling ... ');

    optsPhys.t=optsPhys.tMax;
    optsStocF = optsStocI;
    
    opts.optsPhys = optsPhys;
    opts.optsStoc = optsStocF;

    FCDir = [optsPhys.potNames filesep 'Stochastic' filesep 'Final'];
    
    xFinal = DataStorage(FCDir,@samplepdf,opts,[],~loadSamples);

    fprintf(1,'Finished\n');
    
else
    % return empty data and don't save
    xFinal=[];            
end

%--------------------------------------------------------------------------
% Momentum distribution
%--------------------------------------------------------------------------

if(MBp && doP)

    fprintf(1,'Starting momentum sampling ... ');
    
    optsPhysMB.nParticles = optsPhys.nParticles;
    optsPhysMB.dim        = optsPhys.dim;
    optsPhysMB.m          = optsPhys.m;
    optsPhysMB.kBT        = optsPhys.kBT;
    optsPhysMB.V1DV1      = 'MB';
    optsPhysMB.V2DV2      = 'free';
    optsPhysMB.t          = 0;
    optsPhysMB.geom       = 'full';
    optsPhysMB.type       = 'stoc';
    
    optsStocMB.initialGuess  = 'MBIG';
    optsStocMB.nSamples      = optsStoc.nSamples;
    optsStocMB.thin          = optsStoc.thin;
    optsStocMB.burnin        = optsStoc.burnin;
    
    opts.optsPhys = optsPhysMB;
    opts.optsStoc = optsStocMB;
    
    pDir = [optsPhys.potNames filesep 'Stochastic' filesep 'Momentum'];
    
    pIF = DataStorage(pDir,@samplepdf,opts,[],~loadSamples);

    fprintf(1,'Finished\n');
    
else
    
    fprintf(1,'Using zero momentum ... ');
    
    pIF=zeros(optsStoc.nSamples,optsPhys(1).nParticles*optsPhys(1).dim);
    
    fprintf(1,'Finished\n');
end    


%--------------------------------------------------------------------------
% Dynamics
%--------------------------------------------------------------------------

% placeholder
stocStruct=struct([]);

ICStruct.x = xInitial;
ICStruct.p = pIF;
ICStruct.poolsize = optsStocFull(1).poolsize;

for iStoc=1:nStoc
    % get appropriate options
    optsStoc=optsStocFull(iStoc);

    loadStoc = optsStoc.loadStoc;
    stocName = optsStoc.stocName;
    
    optsStoc = rmfield(optsStoc,'fixedInitial');
    optsStoc = rmfield(optsStoc,'sampleFinal');
    optsStoc = rmfield(optsStoc,'poolsize');
    optsStoc = rmfield(optsStoc,'doStoc');
    optsStoc = rmfield(optsStoc,'loadStoc');
    optsStoc = rmfield(optsStoc,'saveStoc');
    optsStoc = rmfield(optsStoc,'stocName');
    
    opts.optsPhys = optsPhys;
    opts.optsStoc = optsStoc;
     
    fprintf(1,['Starting stochastic dynamics ' num2str(iStoc) '/' num2str(nStoc) ...
                ': ' stocName ' ... ']);

    dynDir = [optsPhys.potNames filesep 'Stochastic' filesep 'Dynamics'];
            
    xpStruct = DataStorage(dynDir,@stochasticStatistics,opts,ICStruct,~loadStoc);
            
    % store in output structure
    stocStruct(iStoc).x=xpStruct.x;
    stocStruct(iStoc).p=xpStruct.p;
    
    stocStruct(iStoc).xInitial=xInitial;
    stocStruct(iStoc).xFinal=xFinal;
    stocStruct(iStoc).pIF=pIF;

    fprintf(1,'Finished\n');
end


