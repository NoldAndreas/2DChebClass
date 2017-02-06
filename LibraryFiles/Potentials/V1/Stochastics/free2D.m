function [VBack_S,VAdd_S]=free2D(x,y,t,optsPhys)

% ensure a column vector (e.g. slicesample uses row vectors)
x=x(:);

VBack=zeros(size(x));

% never use the gradient
DV=zeros(size(x));

% set additional potential to be zero
VAdd=zeros(size(x));


%--------------------------------------------------------------------------

VBack_S = struct('V',VBack, ...
                    'dy1',DV,'dy2',DV, ...
                    'grad', [DV;DV]);

VAdd_S = struct('V',VAdd, ...
                    'dy1',DV,'dy2',DV, ...
                    'grad', [DV;DV]);

end