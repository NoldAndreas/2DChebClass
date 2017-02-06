function [VBack_S,VAdd_S]=bump1D(y,t,optsPhys)
% V = alpha (y-y0)^2

    alpha = optsPhys.alpha;
    beta = optsPhys.beta;
    y0    = optsPhys.y0;

    %--------------------------------------------------------------------------
    VBack  = zeros(size(y));
    DVBack = zeros(size(y));
    %--------------------------------------------------------------------------
    
    if(t==0)
        VAdd  = alpha .* exp(-(y-y0).^2./beta);
        DVAdd = -2*(y-y0)./beta .* VAdd;
    else
        VAdd  = zeros(size(y));
        DVAdd = zeros(size(y));
    end
        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end