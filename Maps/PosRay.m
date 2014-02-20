function [z,dz,dx,ddx,dddx,ddddx] = PosRay(xf,L)
    [z,dz,dx,ddx,dddx,ddddx] = QuotientMap(xf,L,0);        
end
