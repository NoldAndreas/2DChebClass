function [VBack_S,VAdd_S]=oscBump1D(y,t,optsPhys)
% V = alpha0 y^2 + alphaT exp(-(z-z0(t))^2 / beta)

    alpha0 = optsPhys.alpha0;
    alphaT = optsPhys.alphaT;
    beta = optsPhys.beta;
    y0    = 12*sin(t*pi/8);

    %--------------------------------------------------------------------------
    VBack  = alpha0.*y.^2;
    DVBack = 2*alpha0.*y;
    %--------------------------------------------------------------------------
        
    if(t==0)
        VAdd  = zeros(size(y));
        DVAdd = zeros(size(y));
    else
        VAdd  = alphaT .* exp( - (y-y0).^2./beta);
        DVAdd = -2*(y-y0)./beta .* VAdd;
        VAdd(abs(y)==inf) = 0;
        DVAdd(abs(y)==inf) = 0;
    end
        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end