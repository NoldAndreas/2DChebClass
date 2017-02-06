function [meanR,meanFlux,meanV]=getRFluxMeansDDFT(ddft)

rho_t  = ddft.rho_t;
flux_t = ddft.flux_t;

if(isfield(ddft,'v_t'));
    doVelocity = true;
    v_t = ddft.v_t;
end

nTimes=size(rho_t,3);
nSpecies=size(rho_t,2);

% velocity when it exists??
y=ddft.shape.Pts.y;
y(abs(y)==Inf)=0;

Int=ddft.shape.Int;

% mean vectors
meanR    = zeros(nTimes,nSpecies);
meanFlux = zeros(nTimes,nSpecies);  
meanV    = zeros(nTimes,nSpecies);

% calculate mean at each time -- may be able to vectorize this
for iTime=1:nTimes
    % calculate rho and v from full rhov
    rho  = rho_t(:,:,iTime);
    flux = flux_t(:,:,iTime);

    if(doVelocity)
        v = v_t(:,:,iTime);
    end
    
    for iSpecies=1:nSpecies
    
        rhoS  = rho(:,iSpecies);
        fluxS = flux(:,iSpecies);
        
        % calculate (normalised) means
        meanR(iTime,iSpecies)    = ( Int*(rhoS.*y) ) ./ (Int*rhoS);
        meanFlux(iTime,iSpecies) = ( Int*(rhoS.*fluxS) ) ./ (Int*rhoS);
    
        if(doVelocity)
            vS = v(:,iSpecies);
            meanV(iTime,iSpecies) = ( Int*(rhoS.*vS) ) ./ (Int*rhoS);
        end
        
    end
end
    
end