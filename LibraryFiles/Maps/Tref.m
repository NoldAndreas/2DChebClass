function [y,dy,dx,ddx,dddx,ddddx] = Tref(x,d,ep)
%see M1Tref for more capabilities
    c1  = (asinh((1-d)/ep)+asinh((1+d)/ep))/2;
    c2  = asinh((1-d)/ep);
        
    y   = d + ep*sinh(c1*(x-1)+c2);     
    y(x==1)  = 1;
    y(x==-1) = -1;
    
    dy      = ep*c1*cosh(c1*(x-1)+c2);
    
    z       = 1./(ep*(1+((d-y)/ep).^2).^(1/2));
    dx      = z/c1;
    ddx     = (d-y).*z.^3/c1;
    dddx    = 3*((d-y).^2).*(z.^5)/c1 - (z.^3)/c1;
    ddddx   = 15*((d-y).^3).*(z.^7)/c1 - 9*(d-y).*(z.^5)/c1;
    
    

end