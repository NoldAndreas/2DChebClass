function [VBack_S,VAdd_S]=Vext_Exp_HardWall(y1,y2,t,opts)
%Currently designed for HalfSpace

% stoc: inputs are (x,[],t,opts)
% ddft: inputs are (r,theta,t,opts)

    lambda          = 1;
    epsilon_w       = opts.epsilon_w;        
    if(t ~= 0)
        t           = t/opts.tau;
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
    
    alpha      = 2*pi*(1/(2*lambda*exp(lambda))-(1/4)*sqrt(pi)*erf(sqrt(lambda))/lambda^(3/2)+(1/4)*sqrt(pi)/lambda^(3/2));    
    c          = (-16/9*pi)/alpha;    
    	
    y2Wall = y2-0.5 + 1;
	wall   = c*0.5*(pi/lambda)^(3/2)*(1-erf(sqrt(lambda)*y2Wall));       
    dwall  = -c*(pi/lambda)*exp(-lambda*y2Wall.^2);    
    %[wall,dwall] = BarkerHendersonWall(y2 - 0.5,1);
    
    VAdd           = epsilon_w * wall;
    dVAdd          = epsilon_w * dwall;
    
    if(isfield(opts,'epsilon_w_end'))        
        VAdd    = VAdd + wall*(epsilon_w_end - epsilon_w)*(1 -exp(-t^2));
        dVAdd   = dVAdd + dwall*(epsilon_w_end - epsilon_w)*(1 -exp(-t^2));
    elseif(isfield(opts,'epsilon_w_Amplitude'))
         VAdd    = VAdd  + wall*(opts.epsilon_w_Amplitude)*sin(pi*t);
         dVAdd   = dVAdd + dwall*(opts.epsilon_w_Amplitude)*sin(pi*t);
	elseif(isfield(opts,'epsilon_w_max'))        
        if(t<1)
            VAdd    = VAdd  + wall*(opts.epsilon_w_max)*sin(pi*t);
            dVAdd   = dVAdd + dwall*(opts.epsilon_w_max)*sin(pi*t);
        end
    end
    
    if(isfield(opts,'epsilon_w1'))
        w1_steepness = opts.epsilon_w1;
        L1           = opts.L1;
        VAdd    = VAdd + opts.epsilon_w1*(tanh((y1+opts.L1)/opts.w1_steepness)-tanh((y1-opts.L1)/opts.w1_steepness))/2;
        dVAdd   = VAdd - opts.epsilon_w1*(1/2)*((cosh((y+L1)/w1_steepness)).^2-(cosh((-y+L1)/w1_steepness)).^2)./((cosh((y+L1)/w1_steepness)).^2*(cosh((-y+L1)/w1_steepness)).^2*w1_steepness);
    end
    %VAdd    = VAdd + opts.epsilon_w1/4*exp(-(y1/opts.L1).^2);
    VAdd_S      = struct('V',VAdd);
    VAdd_S.dy1  = zeros(size(VAdd));
    VAdd_S.dy2  = dVAdd;

end