function [nR,meanP,xR]=RPhist2D(R,P,nBins,nParticlesS)
%[nR,nP,xR]=RPhist(R,meanP,nBins)
%   determines histogram data for positions and mean momentum in each of
%   the histogram bins
%
% INPUTS: 
%  R            -- (nParticles,dim,nSamples) matrix of positions
%  P            -- (nParticles,dim,nSamples) matrix of momenta
%  nBins        -- number of bins for histogram of form {nBins1,nBins2}
%  nParticlesS  -- number of particles in each species
%
% OUTPUTS:
%  nR            -- (nBins1,nBins2,nSpecies) matrix of the number of particles in each bin
%  meanP         -- (nBins1,nBins2,2,nSpecies) matrix of the mean momentum of each bin
%  xR            -- (nBins1,nBins2,2,nSpecies) matrix of bin centers

nSamples=size(R,3);

nSpecies=length(nParticlesS);

nBins1=nBins(1);
nBins2=nBins(2);

% find start and end points for each species
speciesStart=[1; cumsum(nParticlesS)+1];
speciesEnd=speciesStart(2:end)-1;

nR=zeros(nBins1,nBins2,nSpecies);
meanP=zeros(nBins1,nBins2,2,nSpecies);
xR=zeros(nBins1,nBins2,2,nSpecies);

for iSpecies=1:nSpecies

    RS=R(speciesStart(iSpecies):speciesEnd(iSpecies),:,:);
    PS=P(speciesStart(iSpecies):speciesEnd(iSpecies),:,:);

    RS1=RS(:,1,:);
    RS1=RS1(:);
    RS2=RS(:,2,:);
    RS2=RS2(:);
    
    PS1=PS(:,1,:);
    PS1=PS1(:);
    PS2=PS(:,2,:);
    PS2=PS2(:);
    
    RS1ends=linspace(min(RS1),max(RS1),nBins1+1);
    RS2ends=linspace(min(RS2),max(RS2),nBins2+1);

    RS1h = RS1ends(2)-RS1ends(1);
    RS2h = RS2ends(2)-RS2ends(1);

    RS1mids = (RS1ends(1:end-1)+RS1h/2).';
    RS2mids = (RS2ends(1:end-1)+RS2h/2).';

    RS1boxes = floor( nBins1* ( RS1-min(RS1) ) / (max(RS1)-min(RS1)) ) +1;
    RS1boxes(RS1boxes==nBins1+1)=nBins1;

    RS2boxes = floor( nBins2* ( RS2-min(RS2) ) / (max(RS2)-min(RS2)) ) +1;
    RS2boxes(RS2boxes==nBins2+1)=nBins2;
    
%    RS1boxesA = dsearchn(RS1mids,RS1);
%    RS2boxesA = dsearchn(RS2mids,RS2);
        
    % get number in each box
    nR(:,:,iSpecies)=( accumarray({RS1boxes RS2boxes}, ones(size(RS1)) ) ).';


    % get sum of P's in each box
    meanP(:,:,1,iSpecies)=accumarray({RS1boxes RS2boxes}, PS1);
    meanP(:,:,2,iSpecies)=accumarray({RS1boxes RS2boxes}, PS2);

    nRtemp=nR(:,:,iSpecies);
    nRtemp(nR(:,:,iSpecies)==0)=1;

    meanP(:,:,1,iSpecies)=meanP(:,:,1,iSpecies)./nRtemp;
    meanP(:,:,2,iSpecies)=meanP(:,:,2,iSpecies)./nRtemp;
    
    % normalize to nParticlesS(iSpecies)
    nR(:,:,iSpecies)=nR(:,:,iSpecies)/nSamples/RS1h/RS2h;
    
    
    [xR1,xR2]=meshgrid(RS1mids,RS2mids);
    
    xR(:,:,1,iSpecies)=xR1;
    xR(:,:,2,iSpecies)=xR2;
        
end

end

