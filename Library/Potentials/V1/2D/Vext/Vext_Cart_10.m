function [VBack_S,VAdd_S]=Vext_Cart_10(y1,y2,t,opts)


    epsilon_w       = opts.epsilon_w;        
    if(t ~= 0)
        t           = t/opts.tau; %tau = half a period
    end

    VBack       = zeros(size(y1));

    DVBackDy1   = zeros(size(y1));
    DVBackDy2   = zeros(size(y1));

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);            
    
    if(~isfield(opts,'epsilon_w_end'))
        epsilon_w_end = 2*epsilon_w;
    else
        epsilon_w_end = opts.epsilon_w_end;
    end
    
    
    wall           = pi*(atan(y2)-pi/2+y2./(1+y2.^2));
    wall(y2==inf)  = 0;   
    
    VAdd    = wall*(epsilon_w + (epsilon_w_Amplitude)*sin(pi*t));
    
    
    VAdd_S  = struct('V',VAdd);

end