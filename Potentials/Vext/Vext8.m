function [V,VDiff,VInt] = Vext8(r,theta,intBound,polar)       
        %symmetric in r: Vext(r,theta) = Vext(-r,theta)
        %
        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u        

        V           = 1./((1+r.^2).^2);  
        dVt         = zeros(size(r));%2./((1+r.^2).^2).*cos(2*theta);  
        ddVt        = zeros(size(r));%-4./((1+r.^2).^2).*sin(2*theta);  
        dVr         = -(4*r./((1+r.^2).^3));
        ddVr        = (4*(5*r.^2-1)./((1+r.^2).^4));
        dVrdVt      = zeros(size(r));%-2*(4*r./((1+r.^2).^3)).*cos(2*theta);
        
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
        
    function vi = vextInt(r,polar)
        if(strcmp(polar,'polar'))
            vi = -1./(2*(1+r.^2));
            vi(r == inf) = 0;
        else
            vi          = (r./(1+r.^2) + atan(r))/2.;
            vi(r==inf)  = atan(inf)/2.;
            vi(r==-inf)  = atan(-inf)/2.;
        end
    end

    function vi = vextIntT(t)
        vi = t;
    end


end