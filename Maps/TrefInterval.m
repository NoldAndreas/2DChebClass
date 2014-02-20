function [z,dz,dx,ddx,dddx,ddddx] = TrefInterval(x,y1,y2,y0,ep_y)

    d  = InvLinearMap(y0,y1,y2);
    ep = 2*ep_y/(y2-y1);
    
    L = (y2-y1)/2;
    
    [z1,dz1,dx1,ddx1,dddx1,ddddx1] = Tref(x,d,ep);        
    [z,dzt]  = LinearMap(z1,y1,y2);
    dz       = dz1.*dzt;
    
    dx       = dx1/L;
    ddx      = ddx1/L.^2;
    dddx     = dddx1/L.^3;
    ddddx    = ddddx1/L.^4;
end