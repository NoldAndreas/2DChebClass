function [VBack_S,VAdd_S]=Vext_Pol_2(r,theta,t,optsPhys)

% stoc: inputs are (x,[],t,optsPhys)
% ddft: inputs are (r,theta,t,optsPhys)

    V0   = optsPhys.V0;
    grav = optsPhys.grav;


    x0 = optsPhys.r0;       
        

    VBack   = V0*(r.^2 - 2*x0*r.*cos(theta) + x0^2);
        
    dvdr    = 2*V0*(r-x0*cos(theta));
    dvdt_r  = V0*2*x0*sin(theta);
    
    VBack(r==inf) = inf;
    dvdr(r==inf)  = inf;

    VBack_S = struct('V',VBack,...
                'dy1',dvdr,'dy2',dvdt_r,...
                'grad',[dvdr;dvdt_r]);

    
    xs      = r.*cos(theta);
    VAdd    = grav*xs.*exp(-0.1*r.^2)*(1-exp(-t^2));
    VAdd(r==inf) = 0;
    VAdd_S  = struct('V',VAdd);            
            

end