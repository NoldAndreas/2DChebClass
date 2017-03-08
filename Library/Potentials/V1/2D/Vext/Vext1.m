function [V,VDiff,VInt] = Vext1(r,theta,intBound,polar)       
        %symmetric in r: Vext(r,theta) = Vext(-r,theta)
        %
        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u        

        V           = (r.^2)./(1+r.^2);  
        dVt         = zeros(size(r));  
        dVr         = 2*r./((1+r.^2).^2);
        
        rAsym       = 1;
        drAsym      = 0;
        
        rr          = r(r==inf);
        V(r==inf)   = ones(size(rr))*rAsym;
        dVr(r==inf) = ones(size(rr))*drAsym;
        dVt(r==inf) = zeros(size(rr))*rAsym;
        
        VDiff       = struct('dy1',dVr,'dy2',dVt);
        
        VInt = 0;
        
        if(nargin > 2)            
            VInt = vextInt(intBound.y1_u,polar) - vextInt(intBound.y1_l,polar);
            VInt = VInt *(intBound.y2_u - intBound.y2_l);
        end
        
    function vi = vextInt(r,polar)
        if(strcmp(polar,'polar'))
            vi = (r.^2 - log(1 + r.^2) )/2;
        else
            vi = r - atan(r);
        end
    end
end