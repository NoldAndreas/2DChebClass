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
    
    % find bounds on R data, and add a bit of a buffer to ensure we capture
    % everything
    Rmax=max(RS)+10*eps;
    Rmin=min(RS)-10*eps;

    % calculate (uniform) bin width
    binWidth=(Rmax-Rmin)/(nBins-2);

    % calculate ends of bins
    binEnds=Rmin:binWidth:Rmax; % has length nBins+1
    
    % add buffer bins at the end
    binEnds=[Rmin-binWidth, binEnds, Rmax+binWidth]; %#ok
    nBinsTemp=nBins;

    for iBin=1:nBinsTemp  % loop through bins
        % find which R are in the bin
        binMask=(binEnds(iBin)<=RS & RS < binEnds(iBin+1));
        % count how many
        nR(iBin,iSpecies)=sum(binMask);
        % find mean momentum
        if(nR(iBin,iSpecies)>0)
            meanP(iBin,iSpecies)=sum(PS(binMask))/nR(iBin,iSpecies);
        end
    end
    
    % normalize to nParticlesS(iSpecies)
    nR(:,iSpecies)=nR(:,iSpecies)/nSamples/binWidth;
    
    % find centres of bins
    xRTemp=binEnds+0.5*binWidth;
    xRTemp=xRTemp(1:nBins);  
    xR(:,iSpecies)=xRTemp(:);
    
end
    
end

