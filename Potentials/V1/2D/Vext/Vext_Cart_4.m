function [VBack_S,VAdd_S]=Vext_Cart_4(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)


    v2struct(optsPhys);

    VBack        = -V0*y1;

    DVBackDy1    = -V0*ones(size(y1));
    DVBackDy2    = zeros(size(y1));

VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);

VAdd         = -grav*y1.*(1-exp(-t^2));
VAdd_S       = struct('V',VAdd);

end