function [RoR,RijInv,nij]=makeRoR(x,dim)
% [RoR,RijInv]=makeRoR(x,dim)
%   returns the tensor (r_ij)\otimes(r_ij) / |r_ij|^2 and |r_ij|^(-1)
%
% INPUTS:
%   x           -- (dim*nParticles,1) vector of positions 
%   dim         -- dimension
%
% OUTPUTS:
%   RoR         -- (dim*nParticles,dim*nParticles) tensor
%   RijInv      -- (dim*nParticles,dim*nParticles) tensor


% number of colloid particles
nParticles=length(x)/dim;

% mask to duplicate each row/column dim times
nvec=1:nParticles;             % row vector of column indices
nmask=nvec(ones(dim,1),:);     % dim rows of 1:nParticles (row length)
% nmask is dim copies of  1 2 ... N;

% mask to copy row/columnwise nParticle times
dvec=(1:dim)';                        % column vector of 1:dim
dmask=dvec(:,ones(nParticles,1));     % nParticle columns each 1:dim
% so dmask(:) is 1:dim nParticle times

% -------------------------------------------------------------------------
% first calculate inverse separations
% -------------------------------------------------------------------------

% get separations 
[Rij,nij]=getRij(x,x,dim);     % nParticles x nParticles
Rij=Rij(nmask,nmask);  % dim*nParticles x dim*nParticles by expanding each
                       % element into a dim x dim block
% and inverse separations for non-zero elements
RijInv=zeros(size(Rij));

% this is faster than .^(-1)
RijInv(Rij~=0)=1./Rij(Rij~=0);  % don't invert the diagonal entries

% -------------------------------------------------------------------------
% now calculate RoR
% -------------------------------------------------------------------------

% each column is now a particle
y=reshape(x,dim,nParticles);

% create arrays to determine (x_i-x_j) = X-Y
% X -- nParticles columns each containing x
% Y -- each column i contain nParticle copies of x_i
% so X-Y is analogous to X-X'
X=x(:,ones(nParticles,1));            % dim*nParticles x nParticles
Y=y(dmask,:);                         % dim*nParticles x nParticles
Z=(X-Y);                              % dim*nParticles x nParticles

% A -- each (i,j)  dim x dim block contains dim vertical x_i-x_j
% B -- each (i,j)  dim x dim block contains dim horizontal x_i-x_j
% so the pointwise product gives (x_i-x_j) \otimes (x_i-x_j) in the (i,j)th
% dim x dim block
A=Z(:,nmask);                  % dim row copies of columns

% take the transpose so we're dealing with row vectors for the particles
% now want to duplicate each row dim times
Zt=Z';                         % nParticles x dim*nParticles      
B=-Zt(nmask,:);                % dim*nParticles x dim*nParticles

% simple but much slower version
% A=kron(Z,ones(1,dim));      % dim*nParticles x dim*nParticles
% B=kron(-Z',ones(dim,1));    % dim*nParticles x dim*nParticles

% find  (x_i-x_j) \otimes (x_i-x_j) / ||x_i-x_j||.^2
%RoR=A.*B./(Rij.^2);
%RoR=A.*B.*RijInv.^2;
% probably faster than .^2
RoR=A.*B.*RijInv.*RijInv;