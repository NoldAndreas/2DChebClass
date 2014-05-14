function [VBack_S,VAdd_S]=Vext_Cart_5(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    V0              = optsPhys.V0;
    epsilon_w       = optsPhys.epsilon_w;    
    y10             = optsPhys.y10;
    y20             = optsPhys.y20;     
    t               = t/optsPhys.tau;

    if(V0 == 0)
        VBack       = zeros(size(y1));
        
        DVBackDy1   = zeros(size(y1));
        DVBackDy2   = zeros(size(y1));
    else
        VBack       = V0*((y1-y10).^2 + (y2-y20).^2);

        DVBackDy1   = 2*V0*(y1-y10);
        DVBackDy2   = 2*V0*(y2-y20);
    end

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);            
    
    if(~isfield(optsPhys,'epsilon_w_end'))
        epsilon_w_end = 2*epsilon_w;
    else
        epsilon_w_end = optsPhys.epsilon_w_end;
    end
    
    if(~isfield(optsPhys,'wallAxis') || optsPhys.wallAxis == 1)
        yw = y1;
    else
        yw = y2;
    end
    
    wall           = pi*(atan(yw)-pi/2+yw./(1+yw.^2));
    wall(yw==inf)  = 0;
    wall(yw==-inf) = -pi^2;
    
    VAdd    = wall*(epsilon_w + (epsilon_w_end - epsilon_w)*(1 -exp(-t^2)));
    VAdd_S  = struct('V',VAdd);

end