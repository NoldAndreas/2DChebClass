function r = makeGridPos(optsPhys)
% put particles on a grid of separation 2*sigma.
% used to get initial configuration for sampling

%nParticles=optsPhys.nParticles;

nParticles=sum(optsPhys.nParticlesS);

if(isfield(optsPhys,'sigma'))
    sigma=optsPhys.sigma; 
else
    sigma=1;
end
dim=optsPhys.dim;

% find the length of grid that will hold all the particles
%gridSize=ceil(nParticles^(1/3));
gridSize=ceil(nParticles^(1/dim));

% uniformly spaced grid around zero on which all the particles will fit
x=-gridSize/2:1:gridSize/2;
x=repmat(x,1,dim);

% all possible combinations of points formed from x
r=mycombnk(x,dim);
% narrow down to all unique triples
r=unique(r,'rows');

% distance of each triple from the origin
norm=sum(r.^2,2);

% order points by distance from origin
[~,order]=sort(norm);
r=r(order,:);

% choose the first nParticles of them
r=r(1:nParticles,:);

rMin=min(r);

for iDim=1:dim
    r(:,iDim)=r(:,iDim)-(rMin(iDim)-1/2);
end

% and reshape to give a vector
r=reshape(r',dim*nParticles,1);

sigmaMax=max(sigma);

% scale by twice (arbitrary) the particle size
r=r*2*sigmaMax;
end

