function [VBack_S,VAdd_S]=V1_Rotating_Cart(y1S,y2S,t,optsPhys)

    % get potential parameters
    V0=optsPhys.V0;

    V0r=optsPhys.V0r;
    alphar=optsPhys.alphar;
    rV=optsPhys.rV;

    tau=optsPhys.tau;

    %--------------------------------------------------------------------------

    VBack        = V0.*(y1S.^2 + y2S.^2);

    DVBackDy1    = 2*V0.*y1S;
    DVBackDy2    = 2*V0.*y2S;

    DVBackDy1(abs(y1S)==inf | abs(y2S)==inf) = 0;
    DVBackDy2(abs(y1S)==inf | abs(y2S)==inf) = 0;

    VBack_S = struct('V',VBack,...
            'dy1',DVBackDy1,'dy2',DVBackDy2,...
            'grad',[DVBackDy1;DVBackDy2]);

    %--------------------------------------------------------------------------

    thetaV = mod(t/tau,2*pi);

    VrV    = V0*rV^2;

    xV     = rV*cos(thetaV);
    yV     = rV*sin(thetaV);

    d2             = (y1S-xV).^2+(y2S-yV).^2;

    VAdd           = V0r.*d2.*exp(-d2/alphar^2) - VrV;
    VAdd(abs(y1S)==inf |  abs(y2S)==inf) = 0;

    DVAddDy1    = -2*VAdd.*(y1S-xV) + 2*V0r.*(y1S-xV).*exp(-d2/alphar^2);
    DVAddDy2    = -2*VAdd.*(y2S-yV) + 2*V0r.*(y2S-xV).*exp(-d2/alphar^2);

    DVAddDy1(abs(y1S)==inf | abs(y2S)==inf)=0;
    DVAddDy2(abs(y1S)==inf | abs(y2S)==inf)=0;

    VAdd_S  = struct('V',VAdd, ...
                'dy1',DVAddDy1,'dy2',DVAddDy2,...
                'grad',[DVAddDy1;DVAddDy2]);
end