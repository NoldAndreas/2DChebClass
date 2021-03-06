function [D,sumDivD] = fullWall3D(x,optsPhys)
% D = RP(x,optsPhys)
%   returns the Rotne-Prager diffusion tensor for particles at positions x
%
% INPUTS:
%   x           -- (dim*nParticles,1) vector of positions 
%   optsPhys    -- struct of size [1 1] containing
%                  D0            (diffusion constant)
%                  sigmaH        (hydrodynamic diameter)
%                  dim           (dimension)
%
% OUTPUTS:
%   D          -- (dim*nParticles,dim*nParticles) diffusion matrix

% obtain constants
D0=optsPhys.D0;
sigmaH=optsPhys.sigmaH;
dim=optsPhys.dim;

%[DRP,sumDivDRP]=makew22D(x,sigmaH,dim);

DRP = makew2(x,sigmaH,dim);

Dwall = wallRP3D(x,sigmaH);

% calculate Rotne-Prager form
D=diag(D0)*(eye(length(x)) + DRP + Dwall);

%sumDivD = sumDivDRP + sumDivDwall;

sumDivD = [];

end
