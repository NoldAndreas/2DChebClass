function [y,dy] = ExpMap(x,L,a)%dx,ddx,dddx,ddddx
% z = a + ( L*exp(1/(1-x)-1) - exp(-1/2) )
 
    h    = L/(exp(1)-exp(1/4));
    y    = a + h*(exp(1./(1-x).^2) - exp(1/4) );
    dy   = exp(1./(1-x).^2)./((1-x).^3)*2*h;    

%     h    = L/(exp(1)-exp(1/2));
%     y    = a + h*(exp(1./(1-x)) - exp(1/2) );
%     dy   = exp(1./(1-x))./((1-x).^2)*h;    
        
end