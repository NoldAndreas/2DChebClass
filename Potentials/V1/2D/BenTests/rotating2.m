function [VBack_S,VAdd_S]=rotating2(r,theta,t,optsPhys)


    % get potential parameters
    V0=optsPhys.V0;

    V0r=optsPhys.V0r;
    alphar=optsPhys.alphar;
    rV=optsPhys.rV;

    tau=optsPhys.tau;

    %--------------------------------------------------------------------------

    VBack = V0.*r.^2;

    DVBackDr       = 2*V0.*r;
    DVBackDtheta_r = zeros(size(theta));

    %--------------------------------------------------------------------------

    x      = r.*cos(theta);
    y      = r.*sin(theta);

    thetaV = mod(t/tau,2*pi);

    VrV    = V0*rV^2;

    xV     = rV*cos(thetaV);
    yV     = rV*sin(thetaV);

    d2             = (x-xV).^2+(y-yV).^2;

    VAdd           = V0r.*d2.*exp(-d2/alphar^2) - VrV;
    VAdd(r == inf) = 0;

    dd_dr          = 2*((x-xV).*cos(theta) + (y-yV).*sin(theta));
    dd_dt_r        = 2*(-(x-xV).*sin(theta) + (y-yV).*cos(theta)); %already divided by r

    VAddPrime      = V0r.*(1-d2/alphar).*exp(-d2/alphar);

    DVAddDr        = VAddPrime.*dd_dr;
    DVAddDr(isnan(DVAddDr))=0;
    
    DVAddDtheta_r  = VAddPrime.*dd_dt_r;  %DVAddDtheta    = VAddPrime.*dd_dt;
    %DVAddDtheta_r(r==0)=0;

    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack, ...
                        'dy1',DVBackDr,'dy2',DVBackDtheta_r, ...
                        'grad', [DVBackDr;DVBackDtheta_r]);

    VAdd_S = struct('V',VAdd, ...
                        'dy1',DVAddDr,'dy2',DVAddDtheta_r, ...
                        'grad', [DVAddDr;DVAddDtheta_r]);


end