function [VBack_S,VAdd_S]=zeroQuad1D(y,t,optsPhys)
% V = alpha (y-y0)^2

    alpha = optsPhys.alpha;
    y0    = optsPhys.y0;

    %--------------------------------------------------------------------------
    VBack  = zeros(size(y));
    DVBack = zeros(size(y));
    %--------------------------------------------------------------------------
    
    if(t==0)
        VAdd  = alpha * (y-y0).^2;
        DVAdd = 2*alpha*(y-y0);
    else
        VAdd  = zeros(size(y));
        DVAdd = zeros(size(y));
    end
        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end