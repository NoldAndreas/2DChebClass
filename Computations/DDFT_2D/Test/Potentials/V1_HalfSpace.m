function [VBack_S,VAdd_S]=V1_HalfSpace(y1,y2,t,optsPhys)

    V0 	         = optsPhys.V0;
    V0Add        = optsPhys.V0Add;
    tau          = optsPhys.tau;
    sigma1       = optsPhys.sigma1;
    sigma2       = optsPhys.sigma2;

    VBack        = V0.*(y1.^2+y2.^2);
    DVBackDy1    = 2*V0.*y1;
    DVBackDy2    = 2*V0.*y2;

    DVBackDy1(~isfinite(y1) | ~isfinite(y2)) = 0;
    DVBackDy2(~isfinite(y1) | ~isfinite(y2)) = 0;
    
    VBack_S = struct('V',VBack,...
                    'dy1',DVBackDy1,'dy2',DVBackDy2,...
                    'grad',[DVBackDy1;DVBackDy2]);

	tSwitch = 1-exp(-t^2/tau);
    
    
    VAdd    = V0Add.*exp(-y2.^2./sigma2.^2).*exp(-y1.^2./sigma1.^2)*tSwitch;
    VAdd(~isfinite(y1) | ~isfinite(y2)) = 0;
    VAdd_S  = struct('V',VAdd);

end