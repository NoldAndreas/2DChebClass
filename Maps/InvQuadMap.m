function x = InvQuadMap(y,L05,LD)
%[z,dz,dx,ddx,dddx] = QuadMap(x,hS12,LD)
% z = LD*eps^2*x/(1-x^2+eps^2);
%[-1,1]    -> [-LD,LD]
%sqrt(0.5) -> L05
% L,LD > 0

    eps = sqrt(sqrt(2)/(2*(LD/L05-sqrt(2))));
    L   = (sqrt(2)/L05-2/LD)^(-1);
    
    x           = (sqrt(L^2 + 4*y.^2*(1+eps^2))-L)./(2*y);
    x(y==0)     = 0;
    x(y==inf)   = 1;
    x(y==-inf)  = -1;
    
end