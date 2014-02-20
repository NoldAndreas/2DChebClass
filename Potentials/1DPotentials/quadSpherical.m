function [VBack_S,VAdd_S]=quadSpherical(y,t,optsPhys)
% V = alpha (y-y0)^2

    alpha = optsPhys.alpha;

    %--------------------------------------------------------------------------
    VBack  = alpha*y.^2;
    DVBack = 2*alpha*y;
    %--------------------------------------------------------------------------
    
    VAdd  = zeros(size(y));
    DVAdd = zeros(size(y));

        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end