function xpStruct=stochasticStatistics(opts,ICStruct)
% [x,p]=stochasticStatistics(xInitial,optsPhys,optsStoc)
%   performs stochastic dynamics calculations
%
% INPUTS: 
%  optsPhys -- a structure of size [1 1] containing at least:
%         geom          (geometry; 'spherical'/'planar')
%         dim           (dimension; 1 or 3)
%         D0            (diffusion constant)
%         kBT           (temperature)
%         nParticles    (number of particles)
%         sigma         (particle diameter)
%         m             (particle mass)
%         HI            (hydrodynamic interactions true/false)
%         V1DV1         (string identifying potential and gradient
%                        function)
%         [pot params]  (parameters for potential)
%
% optsStoc -- a structure of size [1 1] containing at least
%         tMax          (max time)
%         tSteps        (number of time steps)
%         nRuns         (number of stochastic runs)
%         sampleShift   (runIC(i) = sample(sampleShift+sampleSkip*i))
%         sampleSkip    (                    "                      )
%         poolsize      (number of parallel processors to use)
%         type          ('rv'/'r')
%         plotTimes     (vector of save/plotting times)
%     
% xInitial -- a matrix of size (nSamples, dim x nParticles) where
%             nSamples >= sampleShift + nRuns*SampleSkip.  
%             Rows contain intial position vectors.
%
% pInitial -- a matrix of size (nSamples, dim x nParticles) where
%             nSamples >= sampleShift + nRuns*SampleSkip.  
%             Rows contain intial momentum vectors.
%
% OUTPUTS:
%  x           -- a matrix of positions of size 
%                 (nRuns, dim x nParticles, length(plotTimes))
%  p           -- a matrix of momenta of size 
%                  (nRuns, dim x nParticles, length(plotTimes))

optsPhys = opts.optsPhys;
optsStoc = opts.optsStoc;

xInitial = ICStruct.x;
pInitial = ICStruct.p;

% get stochastic sampling options
sampleShift=optsStoc.sampleShift;
sampleSkip=optsStoc.sampleSkip;
nRuns=optsStoc.nRuns;
tSteps=optsStoc.tSteps;
% and number of particles
nParticles=optsPhys.nParticles;
% and dimension
dim=optsPhys.dim;

% whether to include noise
noise=optsStoc.noise;

% set initial conditions
x0=zeros(dim*nParticles,nRuns);
p0=x0;
for iRun=1:nRuns
    x0(:,iRun)=xInitial(sampleShift+sampleSkip*iRun,:)';
    p0(:,iRun)=pInitial(sampleShift+sampleSkip*iRun,:)';
end


% initial velocities are all zero
%p0=zeros(dim*nParticles,1);

% get positions of times at which to plot/save data
[plotPosMask,~]=getPlotPos(optsStoc);
nPlots=length(find(plotPosMask));

x=zeros(nRuns,dim*nParticles,nPlots);
p=x;

% get pool size (number of processors)
poolsize=ICStruct.poolsize;

% if there's already a pool open, close it and kill all jobs
if(matlabpool('size')>0)
    fprintf(1,'\n');
    matlabpool('close','force','local')
end

if(poolsize>1)
    fprintf(1,'\n');
    % open pool for parallel computing
    matlabpool('open','local',poolsize);
    oldPath = path;
    addpath('Stochastic',['Stochastic' filesep 'HI'],'Potentials');
end

tic

% note this for loop can be done in parallel
parfor iRun=1:nRuns
    if(noise)
        % create stochastic force for this run
        f=randn(dim*nParticles,tSteps); 
    else
        % set it to be zero if we want to test dynamics without noise - useful
        % for setting times etc)
        f=zeros(dim*nParticles,tSteps);
    end
      
    % do dynamics
    [x(iRun,:,:),p(iRun,:,:)]=stochasticDynamics(f,x0(:,iRun),p0(:,iRun),optsPhys,optsStoc,plotPosMask);
end

toc

if(poolsize>1)
    matlabpool close
    path(oldPath);
end

xpStruct.x = x;
xpStruct.p = p;
