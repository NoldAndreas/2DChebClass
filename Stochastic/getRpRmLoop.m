function [rp,rm] = getRpRmLoop(R,nParticles,nSamples)

%x=x(:);

rp = zeros(nParticles, nParticles, nSamples);
rm = rp;

for iSample = 1:nSamples
    x = R(:,iSample);
    
    X=x(:,ones(length(x),1));

    rp(:,:,iSample) = max(X,X');
    rm(:,:,iSample) = min(X,X');

end

end

