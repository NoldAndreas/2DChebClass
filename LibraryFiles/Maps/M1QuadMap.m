function [y,Diff] = M1QuadMap(x,L05,LD)
%[z,dz,dx,ddx,dddx] = SqrtMap(x,L05,LD)
%[-1,1]    -> [-LD,LD]   and   
%sqrt(0.5) -> hS12
%
% for L = inf this becomes
% z = L05*x/(sqrt(1 - x^2))
    eps = sqrt(sqrt(2)/(2*(LD/L05-sqrt(2))));
    L   = (sqrt(2)/L05-2/LD)^(-1);
    
    y           = L*x./(1-x.^2+eps^2);
    Diff.dydx   = L*(1+x.^2+eps^2)./((1-x.^2+eps^2).^2);
    Diff.dyddx  = 2*L*x.*(3*eps^2+x.^2+3)./(1+eps^2-x.^2).^3;

    %Diff.dydx  = s + (x.^2).*(s.^3)/h2;
    %Diff.dyddx = 3*x.*(s.^3)/h2 + 3*(x.^3).*(s.^5)/h2^2;
     
end