function [VBack_S,VAdd_S]=quadraticHill1D(y,t,optsPhys)

    Vm      = optsPhys.Vm;
    Vp      = optsPhys.Vp;
    V0Add   = optsPhys.V0Add;
    sigma2Add = optsPhys.sigma2Add;
    y0      = optsPhys.y0;
    tau    = optsPhys.tau;
%    tSwitch = optsPhys.tSwitch;

    %--------------------------------------------------------------------------
    %VBack  = Vm.*y.^2.*exp(-2*y);
    %DVBack = 2*Vm.*y.*exp(-2*y); - 2*Vm.*y.^2.*exp(-2*y);
    
    ypMask = y>=0;
    ymMask = y<0;
    VBack = zeros(size(y));
    DVBack = zeros(size(y));
    
    VBack(ypMask)  = Vp.*y(ypMask).^2;
    VBack(ymMask)  = Vm.*y(ymMask).^2;
    DVBack(ypMask)  = 2*Vp.*y(ypMask);
    DVBack(ymMask)  = 2*Vm.*y(ymMask);
    
%     VBack  = Vm.*y.^2;
%     DVBack  = 2*Vm.*y;



    %--------------------------------------------------------------------------

    
    T = exp(-t/tau);
    
    V0Add = T*V0Add;
    
    mu = y0;
    
    VAdd  = -V0Add.*exp(-(y-mu).^2/2/sigma2Add);
    DVAdd = -(y-mu)/sigma2Add.*VAdd;
    
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end