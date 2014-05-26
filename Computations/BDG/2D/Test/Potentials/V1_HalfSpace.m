function [VBack_S,VAdd_S]=V1_HalfSpace(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    V0 	         = optsPhys.V0;
    V0Add        = optsPhys.V0Add;
    sigma1       = optsPhys.sigma1;
    sigma2       = optsPhys.sigma2;
    tau          = optsPhys.tau;
    y10          = optsPhys.y10;
    y20          = optsPhys.y20;

%     VBack        = V0.*(y1.^2+y2.^2);
%     DVBackDy1    = 2*V0.*y1;
%     DVBackDy2    = 2*V0.*y2;

    VBack        = -V0.*exp(-(y1.^2+y2.^2));
    DVBackDy1    = -2*VBack.*y1;
    DVBackDy2    = -2*VBack.*y2;

    
%     VBack        = zeros(size(y1));
%     DVBackDy1    = zeros(size(y1));
%     DVBackDy2    = zeros(size(y1));

    VBack_S = struct('V',VBack,...
                    'dy1',DVBackDy1,'dy2',DVBackDy2,...
                    'grad',[DVBackDy1;DVBackDy2]);

	tSwitch = 1-exp(-t^2/tau);
                
    VAdd         = -tSwitch *V0Add .* exp(-(((y1-y10).^2)./sigma1.^2 + ((y2-y20).^2)./sigma2.^2));
    VAdd_S       = struct('V',VAdd);

end