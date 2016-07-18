function [PK,PV,tempA] = getPressure1D(R,P,nBins,nParticlesS)

% R, P = nParticles x nSamples
nSamples=size(R,2);

nSpecies=length(nParticlesS);

% find start and end points for each species
speciesStart=[1; cumsum(nParticlesS)+1];
speciesEnd=speciesStart(2:end)-1;

nR=zeros(nBins,nSpecies);
sumP=nR;
sumP2=nR;
meanP=nR;
xR=nR;
PK = nR;
PV = nR;

for iSpecies=1:nSpecies

    % extract R, P for this species
    RS=R(speciesStart(iSpecies):speciesEnd(iSpecies),:);
    PS=P(speciesStart(iSpecies):speciesEnd(iSpecies),:);
    
    Rsamples = RS;
    
    % make into column vector
    RS=RS(:);
    PS=PS(:);
       
    % determine spatial bins
    binEnds=linspace(min(RS),max(RS),nBins+1);
    binWidth = binEnds(2)-binEnds(1);
    binMids = (binEnds(1:end-1)+binWidth/2).';
    
    % bin mids for output
    xR(:,iSpecies)=binMids;
    
    % determine which bin each R is in
    RSboxes = floor( nBins* ( RS-min(RS) ) / (max(RS)-min(RS)) ) +1;
    RSboxes(RSboxes==nBins+1)=nBins;

    % number of particles in each bin
    nR(:,iSpecies)=( accumarray(RSboxes, ones(size(RS)) ) ).';
    
    % find P and P^2 in each box
    sumP(:,iSpecies)=accumarray(RSboxes, PS);    
    sumP2(:,iSpecies)=accumarray(RSboxes, PS.^2);    

    % determine mean momentum
    nRtemp=nR(:,iSpecies);
    nRtemp(nR(:,iSpecies)==0)=1;    
    meanP(:,iSpecies)=sumP(:,iSpecies)./nRtemp;
    
    % normalize to nParticlesS(iSpecies)
    nR(:,iSpecies)=nR(:,iSpecies)/nSamples/binWidth;
    sumP(:,iSpecies)=sumP(:,iSpecies)/nSamples/binWidth;
    sumP2(:,iSpecies)=sumP2(:,iSpecies)/nSamples/binWidth;
    
    % do we actually want to take off the mean in the box??
    PK(:,iSpecies) = sumP2(:,iSpecies) - 2*meanP(:,iSpecies).*sumP(:,iSpecies) ...
                    + nR(:,iSpecies).*meanP(:,iSpecies).^2;
    
    % get max and min R for each pair as a diagonal block for each sample
    % (can we do this more efficiently?)
    [RSp,RSm] = getRpRm(RS,nParticlesS(iSpecies),nSamples);
    
    [RSpL,RSmL] = getRpRmLoop(R,nParticlesS(iSpecies),nSamples);
    
    RijS = getRij(RS,RS,1);
    
    [RijSL,ZijSL] = getRijLoop1D(R);
    
    
    RijSL(:)
    RSpL(:)


    % set anything not in the same sample to be infinitely far apart
    sampleMask = kron(eye(nSamples),ones(nParticlesS(iSpecies)));
    RijS(~sampleMask) = inf;
    
    % replicate for each bin
    RSp   = RSp(:,:,ones(nBins,1));
    RSm   = RSm(:,:,ones(nBins,1));
    RijS  = RijS(:,:,ones(nBins,1));
    
    % find right and left ends of bins
    binR = binEnds(2:end);
    binL = binEnds(1:end-1);
   
    % replicate for particles
    BR(1,1,:) = binR(:);
    BL(1,1,:) = binL(:);
    NS = size(RijS,1);
    BR = BR(ones(NS,1),ones(NS,1),:);
    BL = BL(ones(NS,1),ones(NS,1),:);
    
    
    
    tempA = max(0,RSp - BL);
    tempB = max(0,BR - RSm);
    tempC = min(tempA,tempB);
    tempD = min(tempC,binWidth);
    temp = min(tempD,RijS);
    
    
    
    %sampleMask = sampleMask(:,:,ones(nBins,1));
        
    %temp = temp.*sampleMask
    
    PV = squeeze(sum(sum(temp,1),2))/nSamples/binWidth;
    
    
    % still need to include other terms
    
end
    
end