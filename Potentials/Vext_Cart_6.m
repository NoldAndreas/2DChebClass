function [VBack_S,VAdd_S]=Vext_Cart_6(y1,y2,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    V0              = optsPhys.V0;

    if(isfield(optsPhys,'Axis') && strcmp(optsPhys.Axis,'y2'))
        VBack        = V0*(y2/optsPhys.L2).^2;
        DVBackDy1    = zeros(size(y1));
        DVBackDy2    = 2*V0*y2/optsPhys.L2^2;
    elseif(isfield(optsPhys,'Axis') && strcmp(optsPhys.Axis,'y1y2'))
        VBack        = V0*((y1/optsPhys.L1).^2+(y2/optsPhys.L2).^2);
        DVBackDy1    = 2*V0*(y1/optsPhys.L1^2);
        DVBackDy2    = 2*V0*(y2/optsPhys.L2^2);
    elseif(V0==0)
        VBack        = zeros(size(y1));
        DVBackDy1    = zeros(size(y1));
        DVBackDy2    = zeros(size(y1));
    else
        VBack        = V0*((y1/optsPhys.L1).^2);
        DVBackDy1    = 2*V0*y1/(optsPhys.L1^2);
        DVBackDy2    = zeros(size(y1));
    end
    
    if(V0 == 0)
        VBack     = zeros(size(y1));
        DVBackDy1 = zeros(size(y1));
        DVBackDy2 = zeros(size(y1));
    end

    VBack_S = struct('V',VBack,...
                'dy1',DVBackDy1,'dy2',DVBackDy2,...
                'grad',[DVBackDy1;DVBackDy2]);            
    
    
    if(isfield(optsPhys,'grav_start') && isfield(optsPhys,'grav_end'))
        
        if(isfield(optsPhys,'T'))
            T = optsPhys.T;
        else
            T = 1;
        end
        if(t == 0)
            t_evol  = optsPhys.grav_start;
        else
            t_evol  = (optsPhys.grav_start + (optsPhys.grav_end-optsPhys.grav_start)*(1-exp(-t/T)));
        end
        if(isfield(optsPhys,'L1add') && isfield(optsPhys,'L2add'))
            L1add = optsPhys.L1add;
            L2add = optsPhys.L2add;
        else
            L1add = 4;
            L2add = 4;
        end
        
        VAdd    = t_evol*exp(-(y1/L1add).^2).*exp(-(y2/L2add).^2);
    else            
        VAdd    = zeros(size(y1));        
    end
    VAdd_S  = struct('V',VAdd);

end