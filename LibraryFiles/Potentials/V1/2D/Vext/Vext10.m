function [V,VDiff,VInt] = Vext10(r,t,intBound,polar)       
        %symmetric in r: Vext(r,theta) = Vext(-r,theta)
        %
        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u       
        %Integration is not correct!!!
        
        x0 = 1;
        
        x  = r.*cos(t)  - x0; y = r.*sin(t);
        rn          = sqrt(x.^2 + y.^2);
        
        s           = 1;
        
        V           = exp(-rn.^2./(2*s.^2));    

        dVt         = -x0*r.*sin(t).*V/s^2;
        ddVt        = -x0*r.*(x0*r.*(cos(t)).^2+cos(t)*s^2-x0*r).*V/s^4;
        
        dVr         = -(r - x0*cos(t)).*V/s^2;
        ddVr        = (r.^2 - s^2-2*x0*r.*cos(t)+x0*(cos(t)).^2).*V/(s^4);
        dVrdVt      = -x0*sin(t).*(s^2-r.^2+x0*r.*cos(t)).*V/(s^4);
        
        rAsym       = 0;
        drAsym      = 0;
        
        rr          = r(r==inf);
        V(r==inf)   = ones(size(rr))*rAsym;
        dVr(r==inf) = ones(size(rr))*drAsym;
        dVt(r==inf) = zeros(size(rr))*rAsym;
        ddVt(r==inf)= zeros(size(rr))*rAsym;
        
        VDiff       = struct('dy1',dVr,'dy2',dVt,'ddy1',ddVr,'ddy2',ddVt,...
            'dy1dy2',dVrdVt);
        
        VInt = 0;
        
        if(nargin > 2)            
            vi = vextInt(intBound.y1_u,polar) - vextInt(intBound.y1_l,polar);
            vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
            VInt = vi * vt;
        end
        
%Integration only for polar case up to infinity!!!
    function vi = vextInt(r,polar)
        if(strcmp(polar,'polar'))
            vi = -s^2*exp(-r.^2./(2*s.^2));
            vi(r == inf) = 0;
        else
            vi = sqrt(pi/2)*s*erf(r/(sqrt(2)*s));
        end
    end

    function vi = vextIntT(t)
        vi = t;
    end


end


    