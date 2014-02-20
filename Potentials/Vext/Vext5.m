function [V,VDiff,VInt] = Vext5(y1,y2,intBound,polar)       
        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u
        
        k = 4;
        V           = (y1.^2).*(sin(k*y2));  
        dVr         = 2*y1.*sin(k*y2);        
        ddVr        = 2*sin(k*y2);        
        dVt         = k*(y1.^2).*cos(k*y2);  
        ddVt        = -k^2*(y1.^2).*sin(k*y2);  
        dVrdVt      = 2*k*y1.*cos(k*y2);       
        
        VDiff       = struct('dy1',dVr,'dy2',dVt,'ddy1',ddVr,'ddy2',ddVt,...
                            'dy1dy2',dVrdVt);
        
        VInt        = 0;
        if(nargin > 2)
            vr = vextIntR(intBound.y1_u,polar) - vextIntR(intBound.y1_l,polar);
            vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
            VInt = vr*vt;
        end
    function vi = vextIntR(r,polar)
        if(strcmp(polar,'polar'))
            vi = (r.^4)/4;
        else
            vi = r.^3/3;
        end        
    end

    function vi = vextIntT(t)
        vi = cos(k*t)/k;
    end

end