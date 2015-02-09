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
    CutD       = 0.19;
    
    if(f_LL1(CutD) - f_L1(CutD) > 1e-9)
        PrintErrorPos(f_LL1(CutD) - f_L1(CutD),'BarkerHenderson_2D');
    end
    
    markG1     = (r >=1);
    markL1     = ((r < 1) & (r > CutD));
    markLL1    = (r <= CutD);
    
    z = zeros(size(r));
    
    rt         = r(markG1);
    z(markG1)  = -3/2*pi*(1-(21/32)./(rt.^6))./(rt.^5);
    z(markL1)  = f_L1(r(markL1));    
    z(markLL1) = f_LL1(r(markLL1));
    
    z = z*epsilon;
    
    dzdr_r     = 0;    
    alpha  = -16/9*pi*epsilon;       
    
    function z = f_L1(r)
        z = 3*sqrt(1-r.^2)./(160*r.^10).*...
                        (-105-70*r.^2-56*r.^4+112*r.^6+64*r.^8) ...
                 -3*asin(r)./(32*r.^11).*(32*r.^6-21);
    end

    function z = f_LL1(r)
        z = -48/55-24/91*r.^2-2/15*r.^4-15/187*r.^6-105/1976*r.^8;
    end
    
 end