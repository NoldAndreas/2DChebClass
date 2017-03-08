function [VBack_S,VAdd_S]=Vext_Cart_2(y1,y2,t,optsPhys)

    V0   = optsPhys.V0;
    grav = optsPhys.grav;
    
    
    if(~isfield(optsPhys,'y10'))
        optsPhys.y10 = 0;
    end
	if(~isfield(optsPhys,'y20'))
        optsPhys.y20 = 0;
    end
    
    y10 = optsPhys.y10;
    y20 = optsPhys.y20;

    VBack       = V0*((y1-y10).^2 + (y2-y20).^2);
    dvdy1       = 2*V0*(y1 - y10);
    dvdy2       = 2*V0*(y2 - y20);

    VBack_S = struct('V',VBack,...
                'dy1',dvdy1,'dy2',dvdy2,...
                'grad',[dvdy1;dvdy2]);
            
    %y10     = 2*exp(-t^2);  
    if(~isfield(optsPhys,'a'))
        optsPhys.a = [1,0];
    end	
    
    if(~isfield(optsPhys,'b'))
        optsPhys.b = [0,1];
    end
	
    a = optsPhys.a;
    b = optsPhys.b;
    
    
    VAdd    = -grav*(b(1)*y1+b(2)*y2).*exp(-a(1)*y1.^2-a(2)*y2.^2)*(1-exp(-t^2));
    
    VAdd((y1==inf) | (y1==-inf)) = 0;
    
    VAdd_S  = struct('V',VAdd);                   
           
end