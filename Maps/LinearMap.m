function [z,dz,dx,ddx,dddx,ddddx] = LinearMap(xf,a,b)

        if(isscalar(a) && isscalar(b))
            d = (b-a)*ones(size(xf));            
        else
            d = b-a;
        end
        
        % Finite linear mapping from [-1,1] to [a,b]        
        z = a + (xf+1).*d/2;        
        % The obvious linear map for the above
        dz  = d/2;
        % This is dz/dx
        
        if (nargout == 2), return; end
        
        dx   = 2./d;
        ddx  = zeros(size(z));
        dddx = zeros(size(z));
        ddddx = zeros(size(z));
end