function [VBack_S,VAdd_S,alpha]=zeroPotential(r,theta,t,optsPhys)

    %--------------------------------------------------------------------------
    VBack          = zeros(size(r));
    DVBackDr       = zeros(size(r));
    DVBackDtheta_r = zeros(size(r));
    %--------------------------------------------------------------------------
    VAdd           = zeros(size(r));
    DVAddDr        = zeros(size(r));        
    DVAddDtheta_r  = zeros(size(r));    
    %--------------------------------------------------------------------------
    alpha          = 0;

    VBack_S = struct('V',VBack, ...
                        'dy1',DVBackDr,'dy2',DVBackDtheta_r, ...
                        'grad', [DVBackDr;DVBackDtheta_r]);

    VAdd_S = struct('V',VAdd, ...
                        'dy1',DVAddDr,'dy2',DVAddDtheta_r, ...
                        'grad', [DVAddDr;DVAddDtheta_r]);


end