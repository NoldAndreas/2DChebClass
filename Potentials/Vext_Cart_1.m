function [VBack_S,VAdd_S]=Vext_Cart_1(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    V0 	         = optsPhys.V0;
    y10          = optsPhys.y10;
    y20          = optsPhys.y20;
    grav         = optsPhys.grav;

    VBack        = V0*((y1-y10).^2 + (y2-y20).^2);

    DVBackDy1    = 2*V0*(y1-y10);
    DVBackDy2    = 2*V0*(y2-y20);

    VBack_S = struct('V',VBack,...
                    'dy1',DVBackDy1,'dy2',DVBackDy2,...
                    'grad',[DVBackDy1;DVBackDy2]);

    if(grav == 0)
        VAdd         = zeros(size(y1));
    else
        VAdd         = -grav*(y1+y2).*(1-exp(-t^2));
    end
    VAdd_S       = struct('V',VAdd);

end