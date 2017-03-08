function [V,VDiff,VInt] = Vext12(y1,y2,intBound)
    %intBound is a structur with:
    %   - y1_l
    %   - y1_u
    %   - y2_l
    %   - y2_u

    c1 = 0.4;
    c2 = 0.1;
    c3 = 0.1;


    V           = tanh((y1-y2*c1-c3)/c2);  
    z           = (y1-y2*c1-c3)/c2;
    dVdy1       = (1-(tanh(z)).^2)/c2;
    dVdy2       = -(1-(tanh(z)).^2)*c1/c2;
    dVddy1      = -2*tanh(z).*(1-(tanh(z).^2))/c2^2;
    dVddy2      = -2*tanh(z).*(1-(tanh(z)).^2)*(c1/c2)^2;
    dVdy1dy2    = 2*tanh(z).*(1-(tanh(z)).^2)*c1/c2^2;

     VDiff      = struct('dy1',dVdy1,'dy2',dVdy2,'ddy1',dVddy1,'ddy2',dVddy2,'dy1dy2',dVdy1dy2);
     
     if(nargout == 3 && nargin == 3)
        VInt       = vextInt(intBound.y1_u,intBound.y2_u) + ...
                  - vextInt(intBound.y1_u,intBound.y2_l) + ...
                  - vextInt(intBound.y1_l,intBound.y2_u) + ...
                  + vextInt(intBound.y1_l,intBound.y2_l);             
     end
         
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