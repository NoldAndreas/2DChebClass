function [nR,meanP,xR]=RPhist(R,P,nBins,nParticlesS)
%[nR,nP,xR]=RPhist(R,meanP,nBins)
%   determines histogram data for positions and mean momentum in each of
%   the histogram bins
%
% INPUTS: 
%  R            -- (nParticles,nSamples) matrix of 1D positions
%  P            -- (nParticles,nSamples) matrix of 1D momenta
%  nBins        -- number of bins for histogram
%  nParticlesS  -- number of particles in each species
%
% OUTPUTS:
%  nR            -- (nBins,nSpecies) matrix of the number of particles in each bin
%  meanP         -- (nBins,nSpecies) matrix of the mean momentum of each bin
%  xR            -- (nBins,nSpecies) matrix of bin centers

nSamples=size(R,2);

nSpecies=length(nParticlesS);

% find start and end points for each species
speciesStart=[1; cumsum(nParticlesS)+1];
speciesEnd=speciesStart(2:end)-1;

nR=zeros(nBins,nSpecies);
meanP=nR;
xR=nR;

for iSpecies=1:nSpecies

    RS=R(speciesStart(iSpecies):speciesEnd(iSpecies),:);
    PS=P(speciesStart(iSpecies):speciesEnd(iSpecies),:);
    
    RS=RS(:);
    PS=PS(:);
       
    binEnds=linspace(min(RS),max(RS),nBins+1);
    
    binWidth = binEnds(2)-binEnds(1);
    
    binMids = (binEnds(1:end-1)+binWidth/2).';
    
    RSboxes = floor( nBins* ( RS-min(RS) ) / (max(RS)-min(RS)) ) +1;
    RSboxes(RSboxes==nBins+1)=nBins;

    nR(:,iSpecies)=( accumarray(RSboxes, ones(size(RS)) ) ).';
    
    meanP(:,iSpecies)=accumarray(RSboxes, PS);    

    nRtemp=nR(:,iSpecies);
    nRtemp(nR(:,iSpecies)==0)=1;
    
    meanP(:,iSpecies)=meanP(:,iSpecies)./nRtemp;
    
    % normalize to nParticlesS(iSpecies)
    nR(:,iSpecies)=nR(:,iSpecies)/nSamples/binWidth;
    
    xR(:,iSpecies)=binMids;
    
end
    
end

