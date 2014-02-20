function [V,VDiff,VInt] = Vext13(r,th,intBound)
    %intBound is a structur with:
    %   - y1_l
    %   - y1_u
    %   - y2_l
    %   - y2_u    
    c = 1;        
    [x,y] = pol2cart(th,r);

    V           = tanh(x/c);  
    z           = x/c;
    dVdx        = (1-(tanh(z)).^2)/c;
    dVddx       = -2*tanh(z).*(1-(tanh(z).^2))/c^2;    
        
    dVdth       = -dVdx.*r.*sin(th);
    dVdr        = dVdx.*cos(th);
    
    dVddth      = -dVdx.*r.*cos(th)+dVddx.*(r.*sin(th)).^2;
    dVddr       = dVddx.*(cos(th)).^2;
    dVdrdth     = -dVdx.*sin(th)-dVddx.*r.*sin(th).*cos(th);
    
    VDiff      = struct('dy1',dVdr,'dy2',dVdth,'ddy1',dVddr,'ddy2',dVddth,'dy1dy2',dVdrdth,'Lap',dVddx);
     
    VInt       = 0;
     %if(nargout == 3)
     %   VInt       = vextInt(intBound.y1_u,intBound.y2_u) + ...
     %             - vextInt(intBound.y1_u,intBound.y2_l) + ...
     %            - vextInt(intBound.y1_l,intBound.y2_u) + ...
     %             + vextInt(intBound.y1_l,intBound.y2_l);             
     %end
         
%         if(nargin > 2)
%             vr = vextIntR(intBound.y1_u,polar) - vextIntR(intBound.y1_l,polar);
%             vt = vextIntT(intBound.y2_u) - vextIntT(intBound.y2_l);
%             VInt = vr*vt;
%         end
%     function vi = vextIntR(r,polar)
%         if(strcmp(polar,'polar'))
%             vi = (r.^2 - log(1 + r.^2) )/2;
%         else
%             vi = r - atan(r);
%         end        
%     end
% 
     function vi = vextInt(y1,y2)
         Vh           = tanh((y1-y2*c1-c3)/c2);  
         vi = 4*DilogMaple((Vh+1)/2) + 2*log(Vh-1).*log(Vh+1) +...
                                    + (log(Vh+1)).^2 +...
                                    - (log(Vh-1)).^2 +...
                                    - 4*log(2)*log(Vh-1) + ...
                                    + 2*(log(2))^2;
         vi = vi*c2^2/(8*c1);
     end
 
    function y = DilogMaple(x)
        y = dilog(1-x);
    end

end