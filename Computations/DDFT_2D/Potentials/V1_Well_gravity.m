function [VBack_S,VAdd_S]=V1_Well_gravity(y1S,y2S,t,optsPhys)

    % get potential parameters
    V0=optsPhys.V0;

    V0add = optsPhys.V0add;
    sigma1Add = optsPhys.sigma1Add;
    sigma2Add = optsPhys.sigma2Add;
    y10   = optsPhys.y10;
    y20   = optsPhys.y20;
    g     = optsPhys.g;
    gcut  = optsPhys.gcut;
    tau=optsPhys.tau;

    %--------------------------------------------------------------------------

    R             = sqrt(y1S.^2 + y2S.^2);
    DRDy1         = y1S./R;
    DRDy1(y1S==0) = 0;
    DRDy2         = y2S./R;
    DRDy2(y2S==0) = 0;

    [w,Dw] = cutoffWeight(R,5,10);
    DwDy1 = Dw.*DRDy1;
    DwDy2 = Dw.*DRDy2;


    VBack        = V0.*R.^2;

    DVBackDy1    = 2*V0.*y1S;
    DVBackDy2    = 2*V0.*y2S;

    DVBackDy1    = DVBackDy1.*w + VBack.*DwDy1;
    DVBackDy2    = DVBackDy2.*w + VBack.*DwDy2;

    VBack = VBack.*w;
    

    DVBackDy1(abs(y1S)==inf | abs(y2S)==inf) = 0;
    DVBackDy2(abs(y1S)==inf | abs(y2S)==inf) = 0;

    VBack_S = struct('V',VBack,...
            'dy1',DVBackDy1,'dy2',DVBackDy2,...
            'grad',[DVBackDy1;DVBackDy2]);

    %--------------------------------------------------------------------------

    t = t/tau;

    V0add = V0add*exp(-t^2);
    
    VAdd = -V0add.*exp((-(y1S-y10).^2./sigma1Add - (y2S-y20).^2./sigma2Add));
    
    VAdd = VAdd + g*y2S*(1-exp(-t^2)).*exp(-y2S/gcut);
    
    VAdd(abs(y1S)==inf |  abs(y2S)==inf) = 0;
    
    DVAddDy1   = -2*(y1S-y10).*VAdd./sigma1Add;
    DVAddDy2   = -2*(y2S-y20).*VAdd./sigma2Add + exp(-y2S/gcut)*g*(1-exp(-t^2)).*(1-y2S/gcut);
    
    DVAddDy1(abs(y1S)==inf | abs(y2S)==inf)=0;
    DVAddDy2(abs(y1S)==inf | abs(y2S)==inf)=0;

    VAdd_S  = struct('V',VAdd, ...
                'dy1',DVAddDy1,'dy2',DVAddDy2,...
                'grad',[DVAddDy1;DVAddDy2]);
end