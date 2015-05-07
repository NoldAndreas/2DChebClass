function [z,dzdr_r,alpha] = BarkerHenderson_2D(r,parameter) 
%[z,dzdr_r,alpha] = Phi2DLongRange(r,parameter) 
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
    
    if((nargin == 2) && isfield(parameter,'r_cutoff'))
        rc = parameter.r_cutoff;
    else
        rc = inf;
    end

    CutD       = 0.49;    
    %if(f_LL1(CutD) - f_L1(CutD) > 1e-9)
	%  PrintErrorPos(f_LL1(CutD) - f_L1(CutD),'BarkerHenderson_2D');
    %end
	if(isstruct(r))
        r = r.y1_kv;
    end
    
    markG1     = (r >=1);
    markL1     = ((r < 1) & (r > CutD));
    markLL1    = (r <= CutD);
    
    z = zeros(size(r));
    
    rt         = r(markG1);
    z(markG1)  = -3/2*pi*(1-(21/32)./(rt.^6))./(rt.^5);
    z(markL1)  = f_L1(r(markL1));    
    z(markLL1) = f_LL1(r(markLL1));
    dzdr_r     = 0;    
    
    alpha      = 1/2*(-32/9*pi + pi^2*(1/rc^3 - 7/(32*rc^9))); %(16*pi*(1/(3*rc^3) - 1/(9*rc^9)) - 32/9*pi )/2; % 1 - (3/rc^3 - 1/rc^9)/2;    
	c          = epsilon*(-16/9*pi)/alpha;

    z          = c*z;
    dzdr_r     = c*dzdr_r;
    alpha      = c*alpha;
        
    function z = f_L1(r)
        z = 3*sqrt(1-r.^2)./(160*r.^10).*...
                        (-105-70*r.^2-56*r.^4+112*r.^6+64*r.^8) ...
                 -3*asin(r)./(32*r.^11).*(32*r.^6-21);
    end

    function z = f_LL1(r)
        z = -48/55-24/91*r.^2-2/15*r.^4-15/187*r.^6-105/1976*r.^8 ...
            -(3/80)*r.^10-(693/25024)*r.^12-(1287/60800)*r.^14 ...
            -(715/43008)*r.^16-(36465/2732032)*r.^18-(138567/12697600)*r.^20-(29393/3244032)*r.^22-(289731/38010880)*r.^24-(3900225/601358336)*r.^26-(1671525/299892736)*r.^28-(5816907/1203765248)*r.^30-(901620585/213540405248)*r.^32-(12964479/3489660928)*r.^34-(6806351475/2069100494848)*r.^36-(1893496275/646392578048)*r.^38;
    end
    
 end