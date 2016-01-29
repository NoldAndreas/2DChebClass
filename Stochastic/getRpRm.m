function [rp,rm] = getRpRm(x,nParticles,nSamples)

x=x(:);

X=x(:,ones(length(x),1));

sampleMask = ( kron(eye(nSamples),ones(nParticles)) == 1);

Xp = -inf(size(X));
Xp(sampleMask) = X(sampleMask);

Xm = inf(size(X));
Xm(sampleMask) = X(sampleMask);

rp = max(Xp,Xp');
rm = min(Xm,Xm');

rp(~sampleMask) = 0;
rm(~sampleMask) = 0;

end

