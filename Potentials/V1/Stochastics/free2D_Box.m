function [VBack_S,VAdd_S,VGeom_S] = free2D_Box(y1S,y2S,t,optsPhys)

% ensure a column vector (e.g. slicesample uses row vectors)
y1S=y1S(:);

VBack=zeros(size(y1S));

% never use the gradient
DV=zeros(size(y1S));

% set additional potential to be zero
VAdd=zeros(size(y1S));


%--------------------------------------------------------------------------

VBack_S = struct('V',VBack, ...
                    'dy1',DV,'dy2',DV, ...
                    'grad', [DV;DV]);

VAdd_S = struct('V',VAdd, ...
                    'dy1',DV,'dy2',DV, ...
                    'grad', [DV;DV]);

VGeom_S = V1_Box(y1S,y2S,t,optsPhys);
                
end