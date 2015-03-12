function [z,dzdr_r,alpha] = Exponential(r,parameter) 
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

    if(isstruct(r))
        r = r.y1_kv;
    end
    
    markG1     = (r >=1);
    markL1     = (r < 1);    
    
    z          = zeros(size(r));
    
    rt         = r(markG1);
    z(markG1)  = -epsilon*exp(-rt);
    
    dzdr_r     = [];    
    alpha      = -2*pi*exp(-1)*epsilon;
    
 end