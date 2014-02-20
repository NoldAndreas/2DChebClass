function [VBack_S,VAdd_S]=Vext_Pol_4(r,th,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)


V0   = optsPhys.V0;
grav = optsPhys.grav;


erf_v        = (1-erf(r - optsPhys.r0))/2;

%VBack        = -V0*exp(-(r/optsPhys.r0).^2);
VBack        = -V0*erf_v;

DVBackDy1    = 2*V0*r.*exp(-r.^2);
DVBackDy2    = zeros(size(r));

VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);

x = r.*cos(th);

VAdd         = -VBack*(1-exp(-t^2));
VAdd(r == inf) = 0;
VAdd_S       = struct('V',VAdd);

end