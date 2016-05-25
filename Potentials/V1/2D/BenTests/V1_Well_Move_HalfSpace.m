function [VBack_S,VAdd_S,VGeom_S]=V1_Well_Move_HalfSpace(y1S,y2S,t,optsPhys)

    % get potential parameters
    V0=optsPhys.V0;

    V0add = optsPhys.V0add;
    sigma1Add = optsPhys.sigma1Add;
    sigma2Add = optsPhys.sigma2Add;
    y10a   = optsPhys.y10a;
    y20a   = optsPhys.y20a;
    y10b   = optsPhys.y10b;
    y20b   = optsPhys.y20b;
    tau=optsPhys.tau;
    
    n = optsPhys.nAdd;

    %--------------------------------------------------------------------------

    R             = sqrt(y1S.^2 + y2S.^2);
%     DRDy1         = y1S./R;
%     DRDy1(y1S==0) = 0;
%     DRDy2         = y2S./R;
%     DRDy2(y2S==0) = 0;
% 
%     [w,Dw] = cutoffWeight(R,5,10);
%     DwDy1 = Dw.*DRDy1;
%     DwDy2 = Dw.*DRDy2;


    nB = 2;

    VBack       = V0.*R.^nB;
    
    DVBackDy1    = nB*V0*R.^(nB-2).*y1S;
    DVBackDy2    = nB*V0*R.^(nB-2).*y2S;
    

%     DVBackDy1    = DVBackDy1.*w + VBack.*DwDy1;
%     DVBackDy2    = DVBackDy2.*w + VBack.*DwDy2;
% 
%     VBack = VBack.*w;
    

    DVBackDy1(abs(y1S)==inf | abs(y2S)==inf) = 0;
    DVBackDy2(abs(y1S)==inf | abs(y2S)==inf) = 0;

    VBack_S = struct('V',VBack,...
            'dy1',DVBackDy1,'dy2',DVBackDy2,...
            'grad',[DVBackDy1;DVBackDy2]);

    %--------------------------------------------------------------------------
    
    t = t/tau;
    
    %n = 6;
    
    VAdda = - V0add.*exp(-(y1S-y10a).^n./sigma1Add - (y2S-y20a).^n./sigma2Add);
    VAddb = - V0add.*exp(-(y1S-y10b).^n./sigma1Add - (y2S-y20b).^n./sigma2Add);
    
    VAdd = exp(-t^2)*VAdda + (1-exp(-t^2))*VAddb;
    
    VAdd(abs(y1S)==inf |  abs(y2S)==inf) = 0;
    VAdd(abs(y2S)==0) = 0;
    
    DVAddDy1a   = -n*(y1S-y10a).^(n-1).*VAdda./sigma1Add;
    DVAddDy2a   = -n*(y2S-y20a).^(n-1).*VAdda./sigma2Add;
    
    DVAddDy1b   = -n*(y1S-y10b).^(n-1).*VAddb./sigma1Add;
    DVAddDy2b   = -n*(y2S-y20b).^(n-1).*VAddb./sigma2Add;

    DVAddDy1 = exp(-t^2)*DVAddDy1a + (1-exp(-t^2))*DVAddDy1b;
    DVAddDy2 = exp(-t^2)*DVAddDy2a + (1-exp(-t^2))*DVAddDy2b;
    
    DVAddDy1(abs(y1S)==inf | abs(y2S)==inf)=0;
    DVAddDy2(abs(y1S)==inf | abs(y2S)==inf)=0;
    
    VAdd_S  = struct('V',VAdd, ...
                'dy1',DVAddDy1,'dy2',DVAddDy2,...
                'grad',[DVAddDy1;DVAddDy2]);
            
    %----------------------------------------------------------------------

    [VGeom_S] = V1_HardWall(y1S,y2S,t,optsPhys);

end