function [z,dzdr_r,alpha] = Phi1DLongRange(r,parameter) 
    %z      = the 2-particle interaction, distance r
    %dzdr_r = 1/r * dz/dr
    %alpha  = 1/2*2*pi*int( r*f(r), r = 0..infinity )
    
    if(nargin == 1)
        epsilon = 1;
    elseif(isstruct(parameter))
        epsilon = parameter.epsilon;
    else
        epsilon = parameter;
    end

    z      = -2*pi./((1+r.^2).^2)*epsilon;
    dzdr_r = -8*pi*r./((1+r.^2).^3)*epsilon;
    alpha  = -(pi^2/2)*epsilon;
    
 end