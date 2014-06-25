function [z,dzdr_r,alpha] = Phi2DLongRange(r,parameter) 
%[z,dzdr_r,alpha] = Phi2DLongRange(r,parameter) 
%
%z      = -3/2*pi./((1+r.^2).^(5/2))*epsilon;
%dzdr_r = 1/r * dz/dr
%alpha  = -(pi^2/2)*epsilon = 1/2*( 2*pi*int( r*f(r), r = 0..infinity ))
    
    if(nargin == 1)
        epsilon = 1;
    elseif(isstruct(parameter))
        epsilon = parameter.epsilon;
    else
        epsilon = parameter;
    end

    z      = -3/2*pi./((1+r.^2).^(5/2))*epsilon;
    dzdr_r = -3/2*pi*(-5)./((1+r.^2).^(7/2))*epsilon;
    alpha  = -(pi^2/2)*epsilon;
    
 end