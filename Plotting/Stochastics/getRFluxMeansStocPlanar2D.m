function [meanR,meanFlux]=getRFluxMeansStocPlanar2D(x,p,nParticlesS,mS)
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


dim=2;

% get sizes
nRuns=size(x,1);
nParticles=size(x,2)/dim;
nTimes=size(x,3);

% temporary matrix to store 2-dimensional coordinates
R=zeros(nTimes,nParticles,dim,nRuns);
P=R;

for iTime=1:nTimes
    xTemp=squeeze(x(:,:,iTime))';
    pTemp=squeeze(p(:,:,iTime))';

    [RTemp,PTemp]=getRP(xTemp,pTemp,'planar2D',dim);
    
    R(iTime,:,:,:)=RTemp;
    P(iTime,:,:,:)=PTemp;
end
% R,P of size nTimes x nParticles x dim x nRuns

% mean over runs and particles
nSpecies=length(nParticlesS);

% set up output matrices
meanR=zeros(nTimes,nSpecies,dim);
meanFlux=meanR;

% find start and end points for each species
speciesStart=[1; cumsum(nParticlesS)+1];
speciesEnd=speciesStart(2:end)-1;

% calculate mean for each species
for iSpecies=1:nSpecies
    
    RS = R(:,speciesStart(iSpecies):speciesEnd(iSpecies),:,:);
    PS = P(:,speciesStart(iSpecies):speciesEnd(iSpecies),:,:);
    
    % average over runs and then particles
    meanR(:,iSpecies,:)=mean( mean(RS,4), 2);
    meanFlux(:,iSpecies,:)=mean( mean(PS.*RS,4), 2);
    
    meanFlux(:,iSpecies,:)=meanFlux(:,iSpecies,:)/mS(iSpecies);

end