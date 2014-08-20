function [VBack_S,VAdd_S]=Vext_Cart_3(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)


    V0   = optsPhys.V0;
    
    if(isfield(optsPhys,'y10') && isfield(optsPhys,'y20'))
        y10  = optsPhys.y10;
        y20  = optsPhys.y20;
    else
        y10 = 0;
        y20 = 0;
    end
    grav = optsPhys.grav;


    VBack        = V0*((y1-y10).^2 + (y2-y20).^2);

    DVBackDy1    = 2*V0*(y1-y10);
    DVBackDy2    = 2*V0*(y2-y20);

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);

    VAdd         = -grav*y1.*(1-exp(-t^2));
    VAdd_S       = struct('V',VAdd);

end