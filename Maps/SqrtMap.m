function [z,dz,dx,ddx,dddx,ddddx,ddz] = SqrtMap(x,L05,LD)
%[z,dz,dx,ddx,dddx] = SqrtMap(x,L05,LD)
%[-1,1]    -> [-LD,LD]   and   
%sqrt(0.5) -> L05
%
% for L = inf this becomes
% z = L05*x/(sqrt(1 - x^2)) 
    if(nargin == 2)
        LD = inf;
    end        

    h2= L05^2/(  1 - 2*L05^2/LD^2);
    
    z   = x./(sqrt( (1 - x.^2)/h2 + 1/LD^2) );
    s   = ( (1-x.^2)/h2 + 1/LD^2 ).^(-1/2);
    
    dz  = s + (x.^2).*(s.^3)/h2;
    ddz = 3*x.*s.^3/h2 + 3*x.^3.*s.^5/h2^2;
    %z   = hS12*x./((1-x.^2).^(1/2));
    %dz  = hS12./((1-x.^2).^(3/2));
    
    %Inverse
    
    C      = sqrt(1+h2/LD^2);
    %x     = C*z.*(h2+z.^2).^(-1/2);
    dx     = C*h2.*(h2+z.^2).^(-3/2);
    ddx    = -3*C*h2*z.*(h2+z.^2).^(-5/2);
    ddx(z==inf)  = 0;
    ddx(z==-inf) = 0;
    
    dddx   = -3*C*h2*(-4*z.^2+h2).*(h2+z.^2).^(-7/2);        
    dddx(z==inf)  = 0;
    dddx(z==-inf) = 0;    
    
    ddddx  = 15*C*h2*z.*(-4*z.^2+3*h2).*(h2+z.^2).^(-9/2);            
    ddddx(z==inf)  = 0;
    ddddx(z==-inf) = 0;    
     
end