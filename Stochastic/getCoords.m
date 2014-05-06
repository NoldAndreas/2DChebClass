function C=getCoords(x,dim,geom)
%  C=getCetCoords(x,dim,geom)
%   returns a vector of the appropriate d-dimensional position
%   for a vector x in dimension dim and geometry geom (with dimension d).
%
% INPUTS:
%   x           -- (dim*nParticles,1) vector of positions 
%   dim         -- dimension, should be 1, 2 or 3
%   geom        -- geometry
%
% OUTPUTS:
%   C           -- (nParticles,d) vector of d-dimensional positions


% ensure x is a column vector (useful for e.g. slicesample which outputs a
% row vector)
x=x(:);

% calculate number of particles
nParticles=length(x)/dim;

% reshape so each column holds a particle
x=reshape(x,dim,nParticles);

switch geom         % for different geometries

    case 'spherical'   % dim = 3, d = 1
        % R should be the L^2 norm which is now sqrt of sum square of 
        % column elements 
        C=sqrt(sum(x.^2,1));
        C=C(:);
        
    case 'planar'     % dim = 1,2,3, d=1
        % R should be the dim-th coordinate
        C=x(dim,:);
        % reform into column vector
        C=C(:);
        
    case 'planar2D'   % dim = 2, d=2
        % in 2D, it should be a column vector of length(nParticles) for each
        % coordinate.  Then we can use, e.g. x1=C(:,1), y2=C(:,2)
        
        C=x';
            
    case 'polar2D'    % dim = 2, d=2
        % column vector of length(x) for r and theta
        x1=x(1,:);
        x2=x(2,:);
        x1=x1(:);
        x2=x2(:);
        
        [thetar,absr]=cart2pol(x1,x2);
        thetar=mod(thetar,2*pi);
        
        C=[absr,thetar];
        
    case 'full'
        
        C=x(:);
        
end