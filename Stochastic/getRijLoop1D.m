function [Rij,Zij]=getRijLoop1D(R)
% returns unsigned Rij

nParticles = size(R,1);
nSamples   = size(R,2);

Zij = zeros(nParticles, nParticles, nSamples);
Rij = zeros(nParticles, nParticles, nSamples);

for iSample = 1:nSamples

    x = R(:,iSample);
    X=x(:,ones(nParticles,1));            % dim*nParticles x nParticles
    
    Zij(:,:,iSample) = X - X';
    Rij(:,:,iSample) = abs(Zij(:,:,iSample));
end
   
end




