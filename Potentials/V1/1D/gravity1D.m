function [VBack_S,VAdd_S]=gravity1D(y,t,optsPhys)

    g = optsPhys.g;

    %--------------------------------------------------------------------------
    VBack           = zeros(size(y));
    DVBack          = zeros(size(y));        
    %--------------------------------------------------------------------------
    VAdd          = g*y;
    DVAdd         = g*ones(size(y));
    %--------------------------------------------------------------------------

    VBack_S = struct('V',VBack,'DV',DVBack);

    VAdd_S = struct('V',VAdd,'DV',DVAdd);

end