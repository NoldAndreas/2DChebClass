function x = InvQuotientMap(z,L0,a,b)
% h0  = 2*L*dz0/(L-2*dz0)    (if L == inf)  h0  = 2*dz0
% eps = 2*dz0/(L-2*dz0)      (if L == inf)  eps = 0

% x = (-h0 + 2*(z-a)*(1+eps))/(h0 + 2*(z-a))

%    if((nargin <= 3) || (b == -inf) || (b == inf))      
        if(nargin == 2)
            a = 0;
        end
        if(nargin <= 3)
            b = inf(size(L0));
        end
%        eps = 0;
%    else
%        eps = h0./(b-a);
%    end
    if(abs(2*L0/(b-a)-1) < 10^(-6))
        x = -1 + 2*(z - a)/(b-a);
        return;
    end
     
    

    mark  = ((nargin <= 3) | (b == -inf) | (b == inf));
    markN = (mark == 0);
    
    h0        = 2*L0;
    L         = b(markN)-a(markN);
    h0(markN) = h0(markN).*L./(L - 2*L0(markN));
    
    eps        = zeros(size(h0)); 
    eps(markN) = h0(markN)./(b(markN)-a(markN));    
    

    zt  =   bsxfun(@plus,z,-a);

    t1  =   bsxfun(@plus,bsxfun(@times,2*zt,1+eps),-h0);
    t2  =   bsxfun(@plus,h0,2*zt);
    
    x   =   t1./t2;
    
    x(zt==inf)  = 1;
    x(zt==-inf) = 1;
    
end