function [VBack_S,VAdd_S]=V1_InfCapillary(y1,y2,t,optsPhys)

    V0   = optsPhys.V0;
    grav = optsPhys.grav;

    VBack       = V0*(y1.^2 + y2.^2);
    dvdy1       = 2*V0*y1;
    dvdy2       = 2*V0*y2;

    VBack_S = struct('V',VBack,...
                'dy1',dvdy1,'dy2',dvdy2,...
                'grad',[dvdy1;dvdy2]);
            
    %y10     = 2*exp(-t^2);                    
    VAdd    = -grav*y2.*exp(-y1.^2)*(1-exp(-t^2));
    VAdd_S  = struct('V',VAdd);                   
           
end