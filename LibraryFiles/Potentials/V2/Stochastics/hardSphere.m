function [V2,V2prime_r]=hardSphere(Rij,optsPhys)
% [V2,DV2]=hardsphere(x,t,optsPhys)
%   returns the potential and gradient of the hard sphere potential for 
%   particles with separations Rij
%
% INPUTS:
%   Rij         -- matrix of interparticle separations
%   optsPhys    -- struct of size [1 1] containing
%                  dim           (dimension; 1 or 3)
%                  kBT           (temperature)
%                  nParticles    (number of particles)
%                  sigma         (particle diameter)
%
% OUTPUTS:
%   V2           -- scalar, either 0 or NaN
%   DV2          -- (dim*nParticles,1) vector gradient of hard spher 
%                     potential

%-------------------------------------------------------------------------%
% Get parameters                                                          %
%-------------------------------------------------------------------------%

% get constants
sigma=optsPhys.sigma;
kBT=optsPhys.kBT;

%-------------------------------------------------------------------------%
% Calculate V2                                                            %
%-------------------------------------------------------------------------%

% add to self-separation to avoid counting these
Rij=Rij+2*diag(diag(sigma));

%minSep=bsxfun(@plus,sigma',sigma)/2;

minSep=sigma;

sepTest=Rij-minSep;

if(any(any(sepTest<0)))
    V2=NaN;
else
    V2=0;
end

% define a mask depending on the separations for which terms we want to
% keep in DV2
closeMask = Rij<2^(1/24)*minSep;

% include extra 1/r from grad(|x|)=x/|x|
V2prime_r=- kBT*( 48*minSep.^48.*Rij.^(-50) - 24*minSep.^24.*Rij.^(-26) );

% exclude those not sufficiently close, or which are the same particle
V2prime_r(~closeMask)=0;

end