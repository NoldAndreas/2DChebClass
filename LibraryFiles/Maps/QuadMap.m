function [y,dy,dx,ddx,dddx,ddddx] = QuadMap(x,L05,LD)
%[z,dz,dx,ddx,dddx] = QuadMap(x,hS12,LD)
% z = LD*eps^2*x/(1-x^2+eps^2);
%[-1,1]    -> [-LD,LD]
%sqrt(0.5) -> L05
%L05,LD > 0
 
    eps = sqrt(sqrt(2)/(2*(LD/L05-sqrt(2))));
    L   = (sqrt(2)/L05-2/LD)^(-1);
    
    y    = L*x./(1-x.^2+eps^2);
    dy   = L*(1+x.^2+eps^2)./((1-x.^2+eps^2).^2);
    
    z    = sqrt(L^2+4*y.^2*(1+eps^2));
    
    %1st derivative
    dx          = L*(z - L)./(2*z.*y.^2);    
    dx(y==0)    = (1+eps^2)/L;
    dx(y==inf | y == -inf)  = 0;
    
    %2nd derivative
    ddx         = L*(2*y.^2.*(3*L-2*z)*(1+eps^2)+L^2*(L-z))./((z.*y).^3);  
    ddx(y==0)   = 0;
    ddx(y==inf | y == -inf)  = 0;
    
    %3rd derivative    
    dddx        = -3*L*(  (1+eps^2)*(   2*y.^2*L^2.*(5*L-4*z) + 16*y.^4.*(2*L-z)*(1+eps^2)    ) + L^4*(L-z)   )./(z.^5.*y.^4);
    dddx(y==0)  = -6*(1+eps^2)^2/L^3;
    dddx(y==inf | y == -inf)  = 0;
    
    %4th derivative
    IE2         = (1+eps^2);
    ddddx       = 12*L*(...
                    IE2^3*y.^6.*(160*L-64*z)+...
                    L^2*IE2^2*y.^4.*(70*L-48*z)+...
                    L^4*IE2*y.^2.*(14*L-12*z) + ...
                    L^6*(L-z))./(y.^5.*z.^7);
    ddddx(y==0)   = 0;
    ddddx(y==inf | y == -inf)  = 0;

end