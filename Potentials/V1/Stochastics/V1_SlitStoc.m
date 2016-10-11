function [VBack_S,VAdd_S]=V1_SlitStoc(y,t,optsPhys)

    sigma = optsPhys.sigma;
    L = optsPhys.L;
    yMin = -L/2;
    yMax = L/2;
    
    
    %--------------------------------------------------------------------------
    VBack          = zeros(size(y));
    VBack(y<yMin + sigma/2) = inf;
    VBack(y>yMax - sigma/2) = inf;
    
    DVBack         = zeros(size(y));
    DVBack(y<yMin + sigma/2) = inf;
    DVBack(y>yMax - sigma/2) = inf;
    
    %--------------------------------------------------------------------------
    VAdd           = zeros(size(y));
    DVAdd          = zeros(size(y));        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end