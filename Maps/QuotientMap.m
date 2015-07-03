function [z,dz,dx,ddx,dddx,ddddx,ddz] = QuotientMap(x,L0,a,b)
%[z,dz,dx,ddx,dddx,ddddx] = QuotientMap(x,L0,a,b)
%[-1,1] -> [a,b]
%     0 -> a + L0
%L0 - x=0 will be mapped to a+L0. L0 should be of the same dimension as
%a,b. It is necessary that L0 < (b-a)/2!!
%a   - left boundary
%b   - right boundary
%
%in case of mapping to the left: h0<0 and a > b.
%default if no b given: positive ray from a
%default if no a,b given: positive ray from zero
%
%Map: z = a + h0/2 * (1+x)/(1-x+eps)
%We require
%
% (1) z(-1) = a
% (2) z(0)  = a + L0
% (3) z(1)  = b
%
% and obtain (for L = b -a): 
%
% h0  = 2*L*L0/(L-2*L0)    (if L == inf):  h0  = 2*L0
% eps = 2*L0/(L-2*L0)      (if L == inf):  eps = 0
%
% The derivatives are:
% dz/dx = h0/2*(2+eps)/(1-x+eps)^2
%
% dx/dz     = 2*h0*(2+eps)/(h0+2*(z-a))^2
% d^2x/dz^2 = -8*h0*(2+eps)/(h0+2*(z-a))^3
% d^3x/dz^3 = 48*h0*(2+eps)/(h0+2*(z-a))^4

    if(nargin == 2)
        a = zeros(size(L0));
    end
    if(nargin == 3)
        b = inf*ones(size(L0));
    end        
    
    if(isscalar(L0) && (abs(2*L0/(b-a)-1) < 10^(-6)))
        [z,dz,dx,ddx,dddx,ddddx] = LinearMap(x,a,b);
        ddz = zeros(size(z));
        return;
    end

    mark  = ((nargin <= 3) | (b == -inf) | (b == inf));
    markN = (mark == 0);
    
    h0        = 2*L0;
    L         = b(markN)-a(markN);
    h0(markN) = h0(markN).*L./(L - 2*L0(markN));
    
    eps        = zeros(size(h0)); 
    eps(markN) = h0(markN)./(b(markN)-a(markN));
    
    %Compute values of function
    
    t1  = 1./bsxfun(@plus,eps,1-x);
    t2  = bsxfun(@times,h0/2,t1);

    z   = bsxfun(@plus,a,bsxfun(@times,1+x,t2));
    dz  = bsxfun(@times,(2+eps).*h0/2,t1.^2);
    ddz = bsxfun(@times,(2+eps).*h0,t1.^3);
    
    if(~isscalar(L0))
        markL = (abs(2*L0./(b-a)-1) < 10^(-6));
        if(sum(markL) > 0)
            [z(markL),dz(markL)] = LinearMap(x(markL),a(markL),b(markL));        
        end
    end
    
    if(~isscalar(a))
        mark0 = ((b-a) == 0);
        z(mark0)  = a(mark0);
        dz(mark0) = 0;
    end
    
    if (nargout == 2), return; end
    
    zt     = bsxfun(@plus,z,-a);
    
    h02    = bsxfun(@plus,h0,2*zt);
    
    dx     = bsxfun(@times,2*h0.*(2+eps),1./(h02.^2));
    ddx    = bsxfun(@times,-8*h0.*(2+eps),1./(h02.^3));
    dddx   = bsxfun(@times,48*h0.*(2+eps),1./(h02.^4));
    ddddx  = bsxfun(@times,-384*h0.*(2+eps),1./(h02.^5));
    
    if(~isscalar(L0))
        markL = (abs(2*L0./(b-a)-1) < 10^(-6));
        [z(markL),dz(markL),dx(markL),ddx(markL),...
            dddx(markL),ddddx(markL)] = ...
                            LinearMap(x(markL),a(markL),b(markL));
        ddz(markL) = zeros(size(z));
        return;
    end
end