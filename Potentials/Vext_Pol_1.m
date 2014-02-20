function [VBack_S,VAdd_S]=Vext_Pol_1(r,th,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)


V0   = optsPhys.V0;
grav = optsPhys.grav;


VBack        = V0*r.^2;

DVBackDy1    = 2*V0*r;
DVBackDy2    = zeros(size(r));

VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);

x = r.*cos(th);
            
VAdd         = grav*x.*exp(-0.1*r.^2)*(1-exp(-t^2));
VAdd(r == inf) = 0;
VAdd_S       = struct('V',VAdd);

end


