function [V,VDiff,VInt] = Vext3(y1,y2,intBound,polar)       
        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u
        
        if(nargin == 3)
            polar = 'cart';
        end
        
        V           = (y1.^2).*cos(y2);  
        dVr         = 2*y1.*cos(y2);        
        dVt         = -(y1.^2).*sin(y2);  
        
        VDiff       = struct('dy1',dVr,'dy2',dVt);
        
        VInt = 0;
        if(nargin > 2)
            vr = vextIntR(intBound.y1_u,polar) - vextIntR(intBound.y1_l,polar);
            vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
            VInt = vr*vt;
        end
    function vi = vextIntR(r,polar)
        if((nargin == 3) || ~strcmp(polar,'polar'))
            vi = r.^3/3;
        else
            vi = (r.^4)/4;
        end        
    end

    function vi = vextIntT(t)
        vi = sin(t);
    end

end