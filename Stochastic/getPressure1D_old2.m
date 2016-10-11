function [PK,PV,binMids,nR,meanP] = getPressure1D(R,P,nBins,optsPhys)

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
    
    %PK(:,iSpecies)./nRtemp

    
    % NEED TO FIX THIS SO IT WORKS OUT V2DV2 PROPERLY FOR EACH SET OF
    % PARTICLES
    
    %get max and min R and Rij for each pair in each sample
    [RSp,RSm] = getRpRmLoop(R,nParticlesS(iSpecies),nSamples);  
    RijS = getRijLoop1D(R);

    % nParticles x nParticles x nSamples
    
%     RijS = RijS(:);
%     RSp = RSp(:);
%     RSm = RSm(:);

    % find right and left ends of bins
    binR = binEnds(2:end);
    binL = binEnds(1:end-1);

    PV = zeros(nBins,1);

    
        
    for iSample = 1:nSamples

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
            temp  = min(tempD,RijSi);   

            PV(iBin) = PV(iBin) - 1/2 * sum(sum(temp.*dPhidr));
        end

    end
    
    PV = PV/nSamples/binWidth;
    
    
%     for iBin = 1:nBins
%         
%         for iSample = 1:nSamples
% 
%             RSpi = RSp(:,:,iSample);
%             RSmi = RSm(:,:,iSample);
%             RijSi = RijS(:,:,iSample);
%     
%             tempA = max(0,RSpi - binL(iBin));
%             tempB = max(0,binR(iBin) - RSmi);
%             tempC = min(tempA,tempB);
%             tempD = min(tempC,binWidth);
%             temp  = min(tempD,RijSi);   
%             
%             [~,dPhidr_r] = V2DV2(RijSi,optsPhys);  
%             dPhidr = dPhidr_r.*RijSi;
%              
%             PV(iBin) = PV(iBin) - 1/2 * sum(sum(temp.*dPhidr));
%         end
%     end
%     
%     PV = PV/nSamples/binWidth;
        
%     % replicate for each bin
%     RSp   = RSp(:,ones(nBins,1));
%     RSm   = RSm(:,ones(nBins,1));
%     RijS  = RijS(:,ones(nBins,1));

    % replicate for each bin
%     RSp   = RSp(:,:,ones(nBins,1));
%     RSm   = RSm(:,:,ones(nBins,1));
%     RijS  = RijS(:,:,ones(nBins,1));

%     % replicate for samples*particles
%     BR = kron(binR,ones(size(RijS,1),1));
%     BL = kron(binL,ones(size(RijS,1),1));

%     % replicate for each bin
%     RSp   = RSp(:,ones(nBins,1));
%     RSm   = RSm(:,ones(nBins,1));
%     RijS  = RijS(:,ones(nBins,1));


%     tempA = max(0,RSp - BL);
%     tempB = max(0,BR - RSm);
%     tempC = min(tempA,tempB);
%     tempD = min(tempC,binWidth);
%     temp  = min(tempD,RijS);    
%     
    
    
    %dPhir_r = zeros(size(
%     [~,dPhidr_r] = V2DV2(RijS(:),optsPhys);  % need to fix this for
                                                %potential with constants
%     dPhidr = reshape(dPhidr_r,size(RijS)).*RijS;

%     PV = - 1/2 * sum(temp.*dPhidr,1)/nSamples/binWidth;
%     PV = PV(:);

    
end
    
end