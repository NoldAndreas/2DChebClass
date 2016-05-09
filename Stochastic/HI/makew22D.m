function [w2,Dw2]=makew22D(x,sigmaH,dim)
% w2=makew2(x,sigmaH,dim)
%   returns part of the Rotne-Prager diffusion tensor for particles at 
%   positions x.  See Rex & Lowen Eur Phys J E 28 139 (2009)
%
% INPUTS:
%   x           -- (dim*nParticles,1) vector of positions 
%   sigmaH      -- hydrodynamic diameter
%   dim         -- dimension
%
% OUTPUTS:
%   w2          -- (dim*nParticles,dim*nParticles) tensor

% get \hhat r_ij \otimes \hat r_ij and r_ij^{-1}
[RoR,RijInv]=makeRoR(x,dim);


% number of particles
nParticles=length(x)/dim;

% mask to copy row/columnwise nParticles times
dvec=(1:dim)';                        % column vector of 1:dim
dmask=dvec(:,ones(nParticles,1));     % nParticles columns each 1:dim
% so dmask(:) is 1:dim nParticles times

% identity matrix in each dim x dim block, except those on the diagonal
eyes=eye(dim);      % matrix for the dim x dim block
eyes=eyes(dmask,dmask) - eye(dim*nParticles);

% quicker than cubing
RijInv3=RijInv.*RijInv.*RijInv;

% sigmaH is length nParticles*dim

sigmaH=bsxfun(@plus,sigmaH,sigmaH')/2;

sigmaHMax = bsxfun(@max,sigmaH,sigmaH');
sigmaHMin = bsxfun(@min,sigmaH,sigmaH');

lambda = sigmaHMax./sigmaHMin;
alpha = (1+lambda.^2)./((1+lambda).^2);

w2=3/8*sigmaH.*RijInv.*(eyes + RoR) + 1/8*alpha.*sigmaH.^3.*RijInv3.*(eyes-3*RoR);


% divergence term
% ONLY FOR ONE SPECIES

[Rij,nij]=getRij(x,x,dim);
RijInv = 1./Rij;


RijInv(isnan(RijInv) | isinf(RijInv)) = 0;

RijInv2 = RijInv.*RijInv;
RijInv4 = RijInv2.*RijInv2;

d = 2;
sigmaH = sigmaH(1,1);
alpha = alpha(1,1);

Dw2_1 = (3-d) *nij(:,:,1) .*( 3/8*sigmaH.*RijInv2 + 3 * 1/8*alpha.*sigmaH.^3.*RijInv4);
Dw2_2 = (3-d) *nij(:,:,2) .*( 3/8*sigmaH.*RijInv2 + 3 * 1/8*alpha.*sigmaH.^3.*RijInv4);

Dw2_1 = sum(Dw2_1,2);
Dw2_2 = sum(Dw2_2,2);

Dw2 = [Dw2_1'; Dw2_2'];
Dw2 = Dw2(:);



