function [VBack_S,VAdd_S]=quadraticGaussian1D(y,t,optsPhys)

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
    
    mu = y0 - 2*y0*(1-T);
    
    VAdd  = -V0Add.*exp(-(y-mu).^2/2/sigma2Add);
    DVAdd = -(y-mu)/sigma2Add.*VAdd;

%     VAdd = zeros(size(y));
%     DVAdd = zeros(size(y));

    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end