function [nR,meanP,xR]=RPhist2D(R,P,nBins,nParticlesS,range)
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

if(nargin<5)
    getRange=true;
else
    getRange=false;
    x1Min=range(1);
    x1Max=range(2);
    x2Min=range(3);
    x2Max=range(4);
end

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
    
    if(getRange)
    
        x1Min=min(RS1);
        x1Max=max(RS1);
        
        x2Min=min(RS2);
        x2Max=max(RS2);

    end

    RS1ends=linspace(x1Min,x1Max,nBins1+1);
    RS2ends=linspace(x2Min,x2Max,nBins2+1);
    
    RS1h = RS1ends(2)-RS1ends(1);
    RS2h = RS2ends(2)-RS2ends(1);

    RS1mids = (RS1ends(1:end-1)+RS1h/2).';
    RS2mids = (RS2ends(1:end-1)+RS2h/2).';

    RS1boxes = floor( nBins1* ( RS1-x1Min ) / (x1Max-x1Min) ) +1;
    RS1boxes(RS1==x1Max)=nBins1;
    
    RS2boxes = floor( nBins2* ( RS2-x2Min ) / (x2Max-x2Min) ) +1;
    RS2boxes(RS2==x2Max)=nBins2;
      
    mask= (RS1boxes>0) & (RS1boxes<=nBins1) & (RS2boxes>0) & (RS2boxes<=nBins2);
    counter=ones(size(RS1));
    
    %RS1Min=min(RS1boxes(mask));
    RS1Max=max(RS1boxes(mask));
    %RS2Min=min(RS2boxes(mask));
    RS2Max=max(RS2boxes(mask));
    
    % flip dimensions due to meshgrid nonsense
    boxMask1=1:RS2Max;
    boxMask2=1:RS1Max;
        
    % get number in each box
    nR(boxMask1,boxMask2,iSpecies)=( accumarray({RS1boxes(mask) RS2boxes(mask)}, counter(mask) ) ).';


    % get sum of P's in each box
    meanP(boxMask1,boxMask2,1,iSpecies)=( accumarray({RS1boxes(mask) RS2boxes(mask)}, PS1(mask)) ).';
    meanP(boxMask1,boxMask2,2,iSpecies)=( accumarray({RS1boxes(mask) RS2boxes(mask)}, PS2(mask)) ).';

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

