function [VBack_S,VAdd_S]=Vext_Cart_9(y1,y2,t,opts)
%Currently designed for HalfSpace

% stoc: inputs are (x,[],t,opts)
% ddft: inputs are (r,theta,t,opts)
    
    epsilon_w       = opts.epsilon_w;    
    if(t ~= 0)
        t               = t/opts.tau;
    end

    if(~isfield(opts,'V0') || (opts.V0 == 0))
        VBack       = zeros(size(y1));
        
        DVBackDy1   = zeros(size(y1));
        DVBackDy2   = zeros(size(y1));
    else
        V0          = opts.V0;
        y10         = opts.y10;
        y20         = opts.y20;     

        
        VBack       = V0*((y1-y10).^2 + (y2-y20).^2);

        DVBackDy1   = 2*V0*(y1-y10);
        DVBackDy2   = 2*V0*(y2-y20);
    end

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);            
    
    if(~isfield(opts,'epsilon_w_end'))
        epsilon_w_end = 2*epsilon_w;
    else
        epsilon_w_end = opts.epsilon_w_end;
    end
    
    
    wall           = pi*(atan(y2)-pi/2+y2./(1+y2.^2));
    wall(y2==inf)  = 0;    %wall(yw==-inf) = -pi^2;
    
    
    VAdd    = wall*(epsilon_w + (epsilon_w_end - epsilon_w)*(1 -exp(-t^2)));    
    VAdd    = VAdd.*(1+opts.epsilon_w1*(1+tanh(y1/opts.w1_steepness))/2);
    %VAdd    = VAdd + opts.epsilon_w1/4*exp(-(y1/opts.L1).^2);
    VAdd_S  = struct('V',VAdd);

end