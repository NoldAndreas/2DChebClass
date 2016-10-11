function [VBack_S,VAdd_S,VGeom_S]=freeBox3D(x,y,z,t,optsPhys)

% ensure a column vector (e.g. slicesample uses row vectors)
x=x(:);

VBack=zeros(size(x));

% never use the gradient
DV=zeros(size(x));

% set additional potential to be zero
VAdd=zeros(size(x));


%--------------------------------------------------------------------------

VBack_S = struct('V',VBack, ...
                    'dy1',DV,'dy2',DV,'dy3',DV, ...
                    'grad', [DV;DV;DV]);

VAdd_S = struct('V',VAdd, ...
                    'dy1',DV,'dy2',DV,'dy3',DV, ...
                    'grad', [DV;DV;DV]);
%--------------------------------------------------------------------------

[VGeom_S] = V1_Box_3D(x,y,z,t,optsPhys);

end