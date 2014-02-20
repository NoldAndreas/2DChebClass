function [z,dzdr_r,alpha] = BarkerHenderson_2D(r,parameter) 
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
    CutD       = 0.1;
    
    markG1     = (r >=1);
    markL1     = ((r < 1) & (r > CutD));
    markLL1    = (r <= CutD);
    
    z = zeros(size(r));
    
    rt         = r(markG1);
    z(markG1)  = -3/2*pi*(1-(21/32)./(rt.^6))./(rt.^5);
    rt         = r(markL1);
    z(markL1)  = 3*sqrt(1-rt.^2)./(160*rt.^10).*...
                        (-105-70*rt.^2-56*rt.^4+112*rt.^6+64*rt.^8) ...
                 -3*asin(rt)./(32*rt.^11).*(32*rt.^6-21);
    rt         = r(markLL1);
    z(markLL1) = -48/55-24/91*rt.^2-2/15*rt.^4-15/187*rt.^6-105/1976*rt.^8;
    
    z = z*epsilon;
    
    dzdr_r     = 0;    
    alpha  = -16/9*pi*epsilon;
    
 end