function [meanR,meanV]=getRVmeansDDFT(ddft)


rho_t=ddft.rho_t;

nTimes=size(rho_t,3);
nSpecies=size(rho_t,2);

% velocity when it exists??

%v_t=ddft.v_t;

y=ddft.shape.Pts.y;
y(abs(y)==Inf)=0;

Int=ddft.shape.Int;

% mean vectors
meanR=zeros(nTimes,nSpecies);
meanV=zeros(nTimes,nSpecies);  


% calculate mean at each time -- may be able to vectorize this
for iTime=1:nTimes
    % calculate rho and v from full rhov
    rho=rho_t(:,:,iTime);
    %v=v_t(:,:,iTime);

    for iSpecies=1:nSpecies
    
        rhoS=rho(:,iSpecies);
        %vS=v(:,iSpecies);
        
        % calculate (normalised) means
        meanR(iTime,iSpecies)=( Int*(rhoS.*y) ) ./ (Int*rhoS);
        %meanV(iTime,iSpecies)=( Int*(rhoS.*vS) ) ./ (Int*rhoS);
    
    end
end
    
end