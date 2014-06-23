function [meanR,meanFlux]=getRFluxMeansDDFT2D(ddft)

rho_t  = ddft.rho_t;
flux_t = ddft.flux_t;

nTimes=size(rho_t,3);
nSpecies=size(rho_t,2);

y1=ddft.shape.Pts.y1_kv;
y2=ddft.shape.Pts.y2_kv;
y1(abs(y1)==Inf)=0;
y2(abs(y2)==Inf)=0;

Int=ddft.shape.Int;

% mean vectors
meanR    = zeros(nTimes,2,nSpecies);
meanFlux = zeros(nTimes,2,nSpecies);  


% calculate mean at each time -- may be able to vectorize this
for iTime=1:nTimes
    % calculate rho and v from full rhov
    rho  = rho_t(:,:,iTime);
    flux = flux_t(:,:,iTime);

    for iSpecies=1:nSpecies
    
        rhoS  = rho(:,iSpecies);
        fluxS = flux(:,iSpecies);
        fluxS1 = fluxS(1:end/2);
        fluxS2 = fluxS(end/2+1:end);
        
        IntRho = Int*rhoS;
        
        % calculate (normalised) means
        meanR(iTime,1,iSpecies)    = ( Int*(rhoS.*y1) ) ./ IntRho;
        meanR(iTime,2,iSpecies)    = ( Int*(rhoS.*y2) ) ./ IntRho;
        meanFlux(iTime,1,iSpecies) = ( Int*(rhoS.*fluxS1) ) ./ IntRho;
        meanFlux(iTime,2,iSpecies) = ( Int*(rhoS.*fluxS2) ) ./ IntRho;
    
    end
end
    
end