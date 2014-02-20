function [z,dz,dx,ddx,dddx,ddddx] = LinearMap(xf,a,b)
        
        % Finite linear mapping from [-1,1] to [a,b]        
        z = a + (xf+1)*(b-a)/2;        
        % The obvious linear map for the above
        dz  = (b-a)/2*ones(size(z));
        % This is dz/dx
        
        if (nargout == 2), return; end
        
        dx   = 2/(b-a)*ones(size(z));
        ddx  = zeros(size(z));
        dddx = zeros(size(z));
        ddddx = zeros(size(z));
end