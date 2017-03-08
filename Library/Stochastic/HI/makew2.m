function w2=makew2(x,sigmaH,dim)
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
%RijInv3=RijInv.^3;

% calculate w2
%w2=3/8*sigmaH*RijInv .*(eyes + RoR) + 1/16*sigmaH^3*RijInv.^3.*(eyes-3*RoR);

% THIS IS PROBABLY WRONG FOR SPECIES WITH DIFFERENT DIAMETERS

%w2=3/8*diag(sigmaH)*RijInv .*(eyes + RoR) + 1/16*diag(sigmaH.^3)*RijInv3.*(eyes-3*RoR);

% sigmaH is length nParticles*dim

sigmaH=bsxfun(@plus,sigmaH,sigmaH')/2;

sigmaHMax = bsxfun(@max,sigmaH,sigmaH');
sigmaHMin = bsxfun(@min,sigmaH,sigmaH');

lambda = sigmaHMax./sigmaHMin;
alpha = (1+lambda.^2)./((1+lambda).^2);

w2=3/8*sigmaH.*RijInv.*(eyes + RoR) + 1/8*alpha.*sigmaH.^3.*RijInv3.*(eyes-3*RoR);


% OLD VERSION -- slow but understandable
% % number of colloid particles
% N=length(x)/dim;
% 
% % define a mask whose columns are the positions of the 3 coordinates of each
% % colloid particle in r
% mask=reshape(1:dim*N,dim,N);
% 
% w2=zeros(dim*N,dim*N);
% % Z2 is a function of |r_i-r_j| and hence is symmetric.  Also, Z2==0 on 
% % diagonal blocks. Only calculate half of it, then add the transpose.
% for iRow=1:N
%     for iCol=iRow+1:N;
%         % r_i-r_j 
%         R=x(mask(:,iRow))-x(mask(:,iCol));
%         
%         % r_{ij} \otimes r_{ij} / || r_{ij} ||^2
%         normR=norm(R);
%         RoR=R*R'/normR^2;
%        
%         w2(mask(:,iRow),mask(:,iCol)) = ( 3/8*sigmaH/normR * ( eye(dim) + RoR )  ...
%             +1/16*(sigmaH/normR)^3 * (eye(dim) - 3*RoR) );
%     end
% end
% % % add the transpose to complete the calculation -- note there are no
% % % diagonal elements so we don't need to worry about double counting
% w2=w2+w2';

