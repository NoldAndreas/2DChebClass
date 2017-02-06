function [meanR,meanV]=getRVmeansStocPlanar2D(x,p,nParticlesS,mS)
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


% size(x)
% [x,p,percentage]=fixxp(x,p,-10,9);
% percentage

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
meanP=meanR;
meanV=meanR;

% find start and end points for each species
speciesStart=[1; cumsum(nParticlesS)+1];
speciesEnd=speciesStart(2:end)-1;

% calculate mean for each species
for iSpecies=1:nSpecies
    
    % extract relevant coordinates
    Rs=R(:,speciesStart(iSpecies):speciesEnd(iSpecies),:,:);
    Ps=P(:,speciesStart(iSpecies):speciesEnd(iSpecies),:,:);

    % avoid bug in matlab 2014a - put runs into first coordinate 
    Rtemp = permute(Rs,[4,1,2,3]);
    Ptemp = permute(Ps,[4,1,2,3]);

    % average over runs
    meanRtemp = mean(Rtemp,1);
    meanPtemp = mean(Ptemp,1);
    
    % put back in correct order
    meanRs = permute(meanRtemp,[2,3,4,1]);
    meanPs = permute(meanPtemp,[2,3,4,1]);
    
    % average over particles
    meanR(:,iSpecies,:)=mean( meanRs, 2);
    meanP(:,iSpecies,:)=mean( meanPs, 2);  


    % average over runs and then particles
%     meanR(:,iSpecies,:)=mean( mean(R(:,speciesStart(iSpecies):speciesEnd(iSpecies),:,:),4), 2);    
%     meanP(:,iSpecies,:)=mean( mean(P(:,speciesStart(iSpecies):speciesEnd(iSpecies),:,:),4), 2);


    % convert to velocities
    meanV(:,iSpecies,:)=meanP(:,iSpecies,:)/mS(iSpecies);

end
