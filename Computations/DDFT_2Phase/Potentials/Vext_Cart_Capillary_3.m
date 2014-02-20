function [VBack_S,VAdd_S]=Vext_Cart_Capillary_3(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    V0              = optsPhys.V0;
    epsilon_w       = optsPhys.epsilon_w;    
    y10             = optsPhys.y10;
    y20             = optsPhys.y20;     
    t               = t/optsPhys.tau;

    VBack        = zeros(size(y1));
    DVBackDy1    = zeros(size(y1));
    DVBackDy2    = zeros(size(y1));

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);            
    
    if(~isfield(optsPhys,'epsilon_w_end'))
        epsilon_w_end = 2*epsilon_w;
    else
        epsilon_w_end = optsPhys.epsilon_w_end;
    end        
    
    wall     = pi*(atan(y1)-pi/2+y1./(1+y1.^2));
    
    y1_2     = y10 - y1;
    wall_2    = pi*(atan(y1_2)-pi/2+y1_2./(1+y1_2.^2));
    
    VAdd    = wall*(epsilon_w + (epsilon_w_end - epsilon_w)*(1 -exp(-t^2))) + ...
               + wall_2 *(epsilon_w_end + (epsilon_w - epsilon_w_end)*(1 -exp(-t^2)));
    VAdd_S  = struct('V',VAdd);

end