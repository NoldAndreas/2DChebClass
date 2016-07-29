function r = makeGridPosScale3D(optsPhys)
% put particles on a grid in a box of size L1 x L2 x L3.
% used to get initial configuration for sampling

nParticles=sum(optsPhys.nParticlesS);

dim=3;

% find the length of grid that will hold all the particles
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

% shift so all positive
rMin=min(r);
for iDim=1:dim
    r(:,iDim)=r(:,iDim)-(rMin(iDim)-1/2);
end

% scale so all particles fit in the box
fudge = 0.9;
scale1 = fudge*optsPhys.L1/max(r(:,1));
scale2 = fudge*optsPhys.L2/max(r(:,2));
scale3 = fudge*optsPhys.L3/max(r(:,3));

r(:,1) = r(:,1)*scale1;
r(:,2) = r(:,2)*scale2;
r(:,3) = r(:,3)*scale3;

% and reshape to give a vector
r=reshape(r',dim*nParticles,1);

end

