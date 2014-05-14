function [VBack_S,VAdd_S]=Vext_Cart_11(y1,y2,t,opts)
%Currently designed for HalfSpace

% stoc: inputs are (x,[],t,opts)
% ddft: inputs are (r,theta,t,opts)

    
    epsilon_w       = opts.epsilon_w;        
    if(t ~= 0)
        t               = t/opts.tau;
    end

    VBack       = zeros(size(y1));
    DVBackDy1   = zeros(size(y1));
	DVBackDy2   = zeros(size(y1));        

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);            
    
    wall           = pi*(atan(y2)-pi/2+y2./(1+y2.^2));
    wall(y2==inf)  = 0;    
    
    if(isfield(opts,'epsilon_w_end'))        
        VAdd    = wall*(epsilon_w + (epsilon_w_end - epsilon_w)*(1 -exp(-t^2)));
    elseif(isfield(opts,'epsilon_w_Amplitude'))
         VAdd    = wall*(epsilon_w + (opts.epsilon_w_Amplitude)*sin(pi*t));
    end
    
    if(isfield(opts,'epsilon_w1'))
        VAdd    = VAdd + opts.epsilon_w1*(tanh((y1+opts.L1)/opts.w1_steepness)-tanh((y1-opts.L1)/opts.w1_steepness))/2;
    end
    %VAdd    = VAdd + opts.epsilon_w1/4*exp(-(y1/opts.L1).^2);
    VAdd_S  = struct('V',VAdd);

end