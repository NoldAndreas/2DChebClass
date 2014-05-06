function V=getV(x,t,optsPhys)
% V=getV(x,t,optsPhys)
%   returns the gradient of the full potential for particles at x at time t
%
% INPUTS:
%   x           -- (dim*nParticles,1) vector of positions 
%   t           -- time
%   optsPhys    -- struct of size [1 1] containing
%                  V1DV1         (string containing V1 potential function)
%                  V2DV2         (string containing V2 potential function)
%                  [pot params]  (parameters for potential)
%                  dim           (dimension; 1 or 3)
%
% OUTPUTS:
%   V          -- (nParticles,1) vector of potential (both 1-
%                  and 2-body) for particles at postions x at time t

x=x(:);

[V1,~]=getV1DV1(x,[],t,optsPhys);

[V2,~]=getV2DV2(x,optsPhys);

% compute the full V
V=V1+V2;

