function Gamma =RPInv(x,optsPhys)
%Gamma = makeGamma(x,optsPhys)
%   returns the correctly weighted inverse of the Rotne-Prager diffusion 
%   tensor for particles at positions x
%
% INPUTS:
%   x           -- (dim*nParticles,1) vector of positions 
%   optsPhys    -- struct of size [1 1] containing
%                  D0            (diffusion constant)
%                  sigmaH        (hydrodynamic diameter)
%                  dim           (dimension)
%                  kBT           (temperature)
%                  m             (mass)
%
% OUTPUTS:
%   Gamma       -- (dim*nParticles,dim*nParticles) friction matrix


% D0=kBT/m/g.  Assume Z1==0.
% Then D=D0(1+\tilde D) and D*Gamma = kBT/m 1 => Gamma = g(1+\tilde D)^(-1)

% obtain constants
D0=optsPhys.D0;
kBT=optsPhys.kBT;
m=optsPhys.m;
sigmaH=optsPhys.sigmaH;
dim=optsPhys.dim;


D=diag(D0)*(eye(length(x))+makew2(x,sigmaH,dim));

% use that Gamma*D=kBT/m 1
Gamma=kBT*diag(m.^(-1))*inv(D);

% first approximation as it removes having to compute the inverse
% NOT POSITIVE DEFINITE
%Gamma=g*( eye(length(x))-makew2(x,sigmaH) );

end

