function [V,VDiff,VInt] = Vext11(r,theta,intBound,polar)       
        %symmetric in r: Vext(r,theta) = Vext(-r,theta)
        %
        %intBound is a structur with:
        %   - y1_l
        %   - y1_u
        %   - y2_l
        %   - y2_u        
        
        minf            = (r==inf | r==-inf);
        rInf            = r(minf); 

        V               = 1./((1+r.^2).^(5/2)).*sin(2*theta);  
        V(minf)         = zeros(size(rInf));        
        
        dVt             = 2./((1+r.^2).^(5/2)).*cos(2*theta);  
        dVt(minf)       = zeros(size(rInf));
        
        ddVt            = -4./((1+r.^2).^(5/2)).*sin(2*theta);  
        ddVt(minf)      = zeros(size(rInf));
        
        dVr             = -(5*r./((1+r.^2).^(7/2))).*sin(2*theta);
        dVr(minf)       = zeros(size(rInf));
                
        ddVr            = (5*(6*r.^2-1)./((1+r.^2).^(9/2))).*sin(2*theta);
        ddVr(minf)      = zeros(size(rInf));
        
        dVrdVt          = -10*(r./((1+r.^2).^(7/2))).*cos(2*theta);
        dVrdVt(minf)    = zeros(size(rInf));
        
        dddVr           = -(105*r.*(2*r.^2-1)./((1+r.^2).^(11/2))).*sin(2*theta);
        dddVr(minf)     = zeros(size(rInf));

        ddddVr          = (105*(16*r.^4-16*r.^2+1)./((1+r.^2).^(13/2))).*sin(2*theta);
        ddddVr(minf)    = zeros(size(rInf));        

        VDiff       = struct('dy1',dVr,'dy2',dVt,'ddy1',ddVr,'ddy2',ddVt,...
            'dy1dy2',dVrdVt,'dddy1',dddVr,'ddddy1',ddddVr);
        
        VInt = 0;
        
        if(nargin > 2)            
            vi = vextInt(intBound.y1_u,polar) - vextInt(intBound.y1_l,polar);
            vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
            VInt = vi * vt;
        end
        
    function vi = vextInt(r,polar)
        if(strcmp(polar,'polar'))
            vi           = -1./(3*(1+r.^2).^(3/2));
            vi(r == inf) = 0;
        else
            vi                     = r.*(3+2*r.^2)./(3*(1+r.^2).^(3/2));
            vi(r==inf | r == -inf) = 2/3;
        end
    end

    function vi = vextIntT(t)
        vi = -cos(2*t)/2;
    end


end