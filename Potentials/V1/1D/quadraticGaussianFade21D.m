function [VBack_S,VAdd_S]=quadraticGaussianFade21D(y,t,optsPhys)

    Vm      = optsPhys.Vm;
    V0Add   = optsPhys.V0Add;
    sigma2Add = optsPhys.sigma2Add;
    y0      = optsPhys.y0;
    tau    = optsPhys.tau;
%    tSwitch = optsPhys.tSwitch;

    %--------------------------------------------------------------------------
    VBack  = Vm.*y.^2;
    DVBack = 2*Vm.*y;
    %--------------------------------------------------------------------------

    T = exp(-t/tau);
    
    V0Add = T*V0Add;
    
    %mu = y0 - 2*y0*(1-T);
    
    mu = y0;
    
    VAdd1  = -V0Add.*exp(-(y-mu).^2/2/sigma2Add);
    DVAdd1 = -(y-mu)/sigma2Add.*VAdd1;
    
    VAdd2  = -V0Add.*exp(-(y+mu).^2/2/sigma2Add);
    DVAdd2 = -(y+mu)/sigma2Add.*VAdd2;

    VAdd = VAdd1 + VAdd2;
    DVAdd = DVAdd1 + DVAdd2;
    
%     VAdd = zeros(size(y));
%     DVAdd = zeros(size(y));

    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end