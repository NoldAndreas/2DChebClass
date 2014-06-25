function Rij=getRij(x,y,dim)
%  Rij=getRij(x,dim)
%   calculates particle separation matrix of size nParticle x nParticles
%
% INPUTS:
%   x           -- (dim*nParticlesx,1) vector of positions 
%   y           -- (dim*nParticlesy,1) vector of positions 
%   dim         -- dimension, should be 1, 2 or 3
%
% OUTPUTS:
%   Rij         -- nParticlesx x nParticlesy matrix of particle separations

% ensure x is a column vector (useful for e.g. slicesample which outputs a
% row vector)
x=x(:);
y=y(:);

% calculate number of particles
nParticlesx=length(x)/dim;
nParticlesy=length(y)/dim;

% reshape so that each column holds a particle
y=reshape(y,dim,nParticlesy);

% two matrices of size dim*nParticles nParticles, one with same vector
% repeated in each row and the other with the vector repeated in each
% column.  So they're the 'transpose' if the vectors are thought of as
% individual elements

% mask to copy row/columnwise nParticle times
dvec=(1:dim)';                        % column vector of 1:dim
dmask=dvec(:,ones(nParticlesx,1));    % nParticle columns each 1:dim
% so dmask(:) is 1:dim nParticle times
X=x(:,ones(nParticlesy,1));            % dim*nParticles x nParticles
Y=y(dmask,:);                         % dim*nParticles x nParticles

% take the difference of each pair of vectors and reshape so that each
% column holds one x_i-x_j.  The norm is then the square root of the sum of
% the column elements squared.
Z=reshape(X-Y,dim,nParticlesx*nParticlesy);

% calculate the norm of each column and reshape to the correct size
Rij=reshape(sqrt(sum(Z.^2,1)),nParticlesx,nParticlesy);

end



