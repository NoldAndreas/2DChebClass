function [VBack_S,VAdd_S]=MB(p,t,optsPhys)

% ensure a column vector (e.g. slicesample uses row vectors)
p=p(:);

% determine mass for each coordinate
m=optsPhys.m;

% and calculate the potential and its gradient, corresponds to
% exp(- |p|^2/(2 m kBT), where the kBT is included in the sampling

% note |p_i|^2 = \sum_{j=1}^dim p_{i,j}^2 so total KE is just summing over
% all coordinates.

VBack=p.^2/2./m;

% never use the gradient
DV=zeros(size(p));

% set additional potential to be zero
VAdd=zeros(size(p));


%--------------------------------------------------------------------------

VBack_S = struct('V',VBack, ...
                    'dy1',DV,'dy2',DV, ...
                    'grad', DV);

VAdd_S = struct('V',VAdd, ...
                    'dy1',DV,'dy2',DV, ...
                    'grad', DV);


end