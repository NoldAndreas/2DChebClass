function [V,VDiff,VInt] = Vext4(y1,y2,intBound,polar)       
        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u
        
        V           = y1.*cos(y2);  
        dVr         = cos(y2);        
        dVt         = -y1.*sin(y2);  
        
        VDiff       = struct('dy1',dVr,'dy2',dVt);
        
        VInt = 0;
        if(nargin > 2)
            vr = vextIntR(intBound.y1_u,polar) - vextIntR(intBound.y1_l,polar);
            vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
            VInt = vr*vt;
        end
    function vi = vextIntR(r,polar)
        if(strcmp(polar,'polar'))
            vi = (r.^3)/3;
        else
            vi = r.^2/2;
        end        
    end

    function vi = vextIntT(t)
        vi = sin(t);
    end

end