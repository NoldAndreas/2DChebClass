function [V,VDiff,VInt] = Vext2(r,theta,intBound,polar)       
        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u
        
        V           = (sin(theta).^2 + 4*cos(theta).^2).*(r.^2)./(1+r.^2);  
        dVt         = -6*sin(theta).*cos(theta).*(r.^2)./(1+r.^2);  
        dVr         = 2*r./((1+r.^2).^2).*(sin(theta).^2 + 4*cos(theta).^2);
        dVtt        = 6*((sin(theta)).^2 - (cos(theta)).^2 ).*(r.^2)./(1+r.^2);  
        dVrr        = (-2)*(sin(theta).^2 + 4*cos(theta).^2).*(-1+3*r.^2)./((1+r.^2).^3);
        
        minf        = (r==inf)|(r==-inf);
        
        tr           = theta(minf);
        V(minf)      = (sin(tr).^2 + 4*cos(tr).^2);
        dVr(minf)    = 0;
        dVt(minf)    = -6*sin(tr).*cos(tr);
        dVrr(minf)   = 0;        
        dVtt(minf)   = 6*((sin(tr)).^2 - (cos(tr)) );          
                
        VDiff       = struct('dy1',dVr,'dy2',dVt,'ddy1',dVrr,'ddy2',dVtt);
                
        VInt = 0;
        if(nargin > 2)
            vr = vextIntR(intBound.y1_u,polar) - vextIntR(intBound.y1_l,polar);
            vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
            VInt = vr*vt;
        end
    function vi = vextIntR(r,polar)
        if(strcmp(polar,'polar'))
            vi = (r.^2 - log(1 + r.^2) )/2;
        else
            vi = r - atan(r);
        end        
    end

    function vi = vextIntT(t)
        vi = 3/2*sin(t)*cos(t) + 5/2*t;
    end

end