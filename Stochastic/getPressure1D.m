function [PK,PV,binMids,nR,meanP,PbyBin] = getPressure1D(R,P,nBins,optsPhys)

nParticlesS = optsPhys.nParticlesS;

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

V2DV2 = str2func(optsPhys.V2DV2);

for iSpecies=1:nSpecies

    % extract R, P for this species
    RS=R(speciesStart(iSpecies):speciesEnd(iSpecies),:);
    PS=P(speciesStart(iSpecies):speciesEnd(iSpecies),:);
    
    % make into column vector
    RS=RS(:);
    PS=PS(:);

    if(~isfield(optsPhys,'rMax'))
        rMax = max(RS);
        rMin = min(RS);
    else
        rMax = optsPhys.rMax;
        rMin = optsPhys.rMin;
    end

    % determine spatial bins
    binEnds=linspace(rMin,rMax,nBins+1);
    binWidth = binEnds(2)-binEnds(1);
    binMids = (binEnds(1:end-1)+binWidth/2).';
    
    % bin mids for output
    xR(:,iSpecies)=binMids;
    
    % add fake particle to each bin to ensure it's counted
    RSTemp = [RS;binMids];
    PSTemp = [PS;zeros(size(binMids))];
    
    % determine which bin each R is in
    RSboxes = floor( nBins* ( RSTemp-rMin ) / (rMax-rMin) ) +1;
        
    % end bins contain all particles outside range so may well be junk
    RSboxes(RSboxes > nBins)=nBins;
    RSboxes(RSboxes < 1)=1;
    
    tempCount = ( accumarray(RSboxes, ones(size(RSTemp)) ) ).';
    tempCount = tempCount - ones(size(tempCount));
    
    % number of particles in each bin
    nR(:,iSpecies) = tempCount;
    
    % find P and P^2 in each box
    sumP(:,iSpecies)=accumarray(RSboxes, PSTemp);    
    sumP2(:,iSpecies)=accumarray(RSboxes, PSTemp.^2);    

    % remove fake particles
    boxMask = RSboxes(1:length(RS));
    
    % determine collection of momenta in each bin
    PbyBin = cell(nBins,1);
    for iBin = 1:nBins
        PbyBin{iBin} = P(boxMask == iBin);
    end
    
    % determine mean momentum
    nRtemp=nR(:,iSpecies);
    nRtemp(nR(:,iSpecies)==0)=1;    
    meanP(:,iSpecies)=sumP(:,iSpecies)./nRtemp;
    
    % normalize to nParticlesS(iSpecies)
    nR(:,iSpecies)=nR(:,iSpecies);
    sumP(:,iSpecies)=sumP(:,iSpecies);
    sumP2(:,iSpecies)=sumP2(:,iSpecies);
    
    PK(:,iSpecies) = sumP2(:,iSpecies) - 2*meanP(:,iSpecies).*sumP(:,iSpecies) ...
                    + nR(:,iSpecies).*meanP(:,iSpecies).^2;
    
    PK(:,iSpecies) = PK(:,iSpecies)/nSamples/binWidth;
    
    nR(:,iSpecies) = nR(:,iSpecies)/nSamples/binWidth;
    
    %get max and min R and Rij for each pair in each sample
    [RSp,RSm] = getRpRmLoop(R,nParticlesS(iSpecies),nSamples);  
    RijS = getRijLoop1D(R);

    % nParticles x nParticles x nSamples

    % find right and left ends of bins
    binR = binEnds(2:end);
    binL = binEnds(1:end-1);

    PV = zeros(nBins,1);

    hw = waitbar(0,'Computing PV');

    for iSample = 1:nSamples

        waitbar(iSample/nSamples,hw);

        RSpi = RSp(:,:,iSample);
        RSmi = RSm(:,:,iSample);
        RijSi = RijS(:,:,iSample);

        [~,dPhidr_r] = V2DV2(RijSi,optsPhys);  
        dPhidr = dPhidr_r.*RijSi;

        for iBin = 1:nBins

            tempA = max(0,RSpi - binL(iBin));
            tempB = max(0,binR(iBin) - RSmi);
            tempC = min(tempA,tempB);
            tempD = min(tempC,binWidth);
            lij  = min(tempD,RijSi);   

            PV(iBin) = PV(iBin) - 1/2 * sum(sum(lij.*dPhidr));
        end

    end

    close(hw);

    PV = PV/nSamples/binWidth;
    
    
    
end
    
end