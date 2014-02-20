function [VBack_S,VAdd_S]=Vext_Pol_2Species_1(r_S,t_S,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    %y1 = r_S(:,1).*cos(t_S(:,1));
    %y2 = r_S(:,1).*sin(t_S(:,1));

    v2struct(optsPhys);

    VBack        = V0*r_S.^2;

    DVBackDy1    = 2*V0*r_S;
    DVBackDy2    = zeros(size(r_S));

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);
            
    t = t/tau;
    %1st Species
    %VAdd(:,1) = (1 -exp(-t^2))*optsPhys.grav1*y1;
    %VAdd(:,2) = (1 -exp(-t^2))*optsPhys.grav2*y1;
    
    y1   = r_S.*cos(t_S);
    if(t == 0)
        VAdd = zeros(size(y1));
    else
        VAdd = (1 -exp(-t^2))*exp(-(r_S/4).^2).*optsPhys.grav.*y1;    
    end
    
    
    VAdd_S  = struct('V',VAdd);

end