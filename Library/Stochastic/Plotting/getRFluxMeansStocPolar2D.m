function [meanR,meanFlux]=getRFluxMeansStocPolar2D(x,p,nParticlesS,mS)

[meanRC,meanFC]=getRFluxMeansStocPlanar2D(x,p,nParticlesS,mS);

dim=2;

nSpecies=length(nParticlesS);
nTimes=size(x,3);

meanR=zeros(nTimes,nSpecies,dim);
%meanV=zeros(nTimes,nSpecies,dim);

for iSpecies=1:nSpecies
    
    meanXS = meanRC(:,iSpecies,1);
    meanYS = meanRC(:,iSpecies,2);

    [meanThetaS,meanRS] = cart2pol(meanXS,meanYS);
    meanThetaS=mod(meanThetaS,2*pi);

    meanR(:,iSpecies,1) = meanRS;
    meanR(:,iSpecies,2) = meanThetaS;
    
end

% NEED TO THINK ABOUT THIS

meanFlux=meanFC;

end