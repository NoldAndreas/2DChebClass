function meanStruct=getMeansParticles(stocStruct)
%meanStruct=getMeansParticles(stocStruct)
% gets mean position and momenta for each particle over stochastic runs.
% Should only be used for distinguishable particles.
%
% INPUTS:
%  stocStruct  --  structure of size [1 1] containing
%                    x       (postiton matrix of size 
%                                     nRun x (dim*nParticles) x nTimes)
%                    p       (momentum matrix of size 
%                                     nRun x (dim*nParticles) x nTimes)
%
% OUTPUTS:
%  meanStruct  -- structure of size [1 1] containing
%                    xMean   (mean position matrix of size
%                                     (dim*nParticles) x nTimes)
%                    pMean   (mean momentum matrix of size
%                                     (dim*nParticles) x nTimes)

% get dimensions
nStoc=size(stocStruct,2);
tempx=stocStruct(1).x;
nRuns=size(tempx,1);

% initialize output structure
meanStruct=struct([]);

for iStoc=1:nStoc
    % get the x and p data for each stochastic calculation
    x=stocStruct(iStoc).x;
    p=stocStruct(iStoc).p;
    
    % calculate mean over runs for each particle and each time
    meanStruct(iStoc).xMean=squeeze(sum(x,1))/nRuns;
    meanStruct(iStoc).pMean=squeeze(sum(p,1))/nRuns;
end
