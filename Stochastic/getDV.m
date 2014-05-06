function DV=getDV(x,t,optsPhys)
% DV=getDV(x,t,optsPhys)
%   returns the gradient of the full potential for particles at x at time t
%
% INPUTS:
%   x           -- (dim*nParticles,1) vector of positions 
%   t           -- time
%   optsPhys    -- struct of size [1 1] containing
%                  V1DV1         (string containing V1 potential function)
%                  [pot params]  (parameters for potential)
%                  dim           (dimension; 1 or 3)
%                  kBT           (temperature)
%                  nParticles    (number of particles)
%                  sigma         (particle diameter)
%
% OUTPUTS:
%   DV          -- (dim*nParticles,1) vector gradient of potential (both 1-
%                  and 2-body) for particles at postions x at time t

% ensure column vector
x=x(:);

% evaluate it for DV1, note the first output is V1
[~,DV1]=getV1DV1(x,[],t,optsPhys);

% evaluate DV2
[~,DV2]=getV2DV2(x,optsPhys);

% compute the full DV
DV=DV1+DV2;