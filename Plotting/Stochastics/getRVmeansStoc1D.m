function [meanR,meanV]=getRVmeansStoc1D(x,p,geom,dim,nParticlesS,mS)
% [meanR,meanV]=getRVmeansStoc(x,p,geom,dim)
%   returns mean position and velocity for stochastic results at each time
%
% INPUTS:
%   x           -- (nSamples,dim*nParticles,nTimes) matrix of positions 
%   p           -- (nSamples,dim*nParticles,nTimes) matrix of momenta 
%   dim         -- dimension, should be 1 or 3
%   geom        -- geometry, should be 'planar' or 'spherical'
%
% OUTPUTS:
%   meanR       -- (nTimes,nSpecies) matrix of mean positions
%   meanP       -- (nTimes,nSpecies) matrix of mean velocities


% get sizes
nRuns=size(x,1);
nParticles=size(x,2)/dim;
nTimes=size(x,3);

% temporary matrix to store one-dimensional coordinates
R=zeros(nRuns,nParticles,nTimes);
P=R;

for iTime=1:nTimes
    xTemp=squeeze(x(:,:,iTime))';
    pTemp=squeeze(p(:,:,iTime))';

    [RTemp,PTemp]=getRP(xTemp,pTemp,geom,dim);
     
    R(:,:,iTime)=RTemp';
    P(:,:,iTime)=PTemp';
end
% R,P of size nRuns x nParticles x nTimes

% mean over runs and particles
nSpecies=length(nParticlesS);

% set up output matrices
meanR=zeros(nTimes,nSpecies);
meanP=meanR;
meanV=meanR;

% find start and end points for each species
speciesStart=[1; cumsum(nParticlesS)+1];
speciesEnd=speciesStart(2:end)-1;

% calculate mean for each species
for iSpecies=1:nSpecies
    
    % extract relevant coordinates
    Rs=R(:,speciesStart(iSpecies):speciesEnd(iSpecies),:);
    Ps=P(:,speciesStart(iSpecies):speciesEnd(iSpecies),:);
    
    % average over runs and then particles
    meanR(:,iSpecies)=squeeze(mean( mean(Rs,1), 2));
    meanP(:,iSpecies)=squeeze(mean( mean(Ps,1), 2));
    % convert to velocities
    meanV(:,iSpecies)=meanP(:,iSpecies)/mS(iSpecies);

end