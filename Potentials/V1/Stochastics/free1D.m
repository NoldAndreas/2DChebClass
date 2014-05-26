function [VBack_S,VAdd_S]=zeroPotential1D(y,t,optsPhys)

    %--------------------------------------------------------------------------
    VBack          = zeros(size(y));
    DVBack         = zeros(size(y));
    %--------------------------------------------------------------------------
    VAdd           = zeros(size(y));
    DVAdd          = zeros(size(y));        
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end