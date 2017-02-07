function [z,dzdr_r,alpha] = Phi2DLongRange(r,parameter) 
%
%z      = -3/2*pi./((1+r.^2).^(5/2))*epsilon;
%dzdr_r = 1/r * dz/dr
%alpha  = -(pi^2/2)*epsilon = 1/2*( 2*pi*int( r*f(r), r = 0..infinity ))

    if(nargin < 2)
        epsilon = 1;
    elseif(isstruct(parameter))
        epsilon = parameter.epsilon;
    else
        epsilon = parameter;
    end

    if(isstruct(r))
        r = r.y1_kv;
    end   

    z      = epsilon*(-3/2)*pi./((1+r.^2).^(5/2));
    dzdr_r = epsilon*(-3/2)*pi*(-5)./((1+r.^2).^(7/2));
    alpha  = -(pi^2/2)*epsilon;
    
 end