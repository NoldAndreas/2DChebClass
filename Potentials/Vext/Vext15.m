function [V,VDiff,VInt] = Vext15(x,y)

    disp('Integration only valid for half-space y2>0');

    a = 0.5;

    V            = exp(-x.^2-a*y.^2);        
    dVdx         = -2*x.*V;
    dVdx(V==0)   = 0;
    dVdy         = -2*a*y.*V;
    dVdy(V==0)   = 0;
    
    dVddx        = (-2+(2*x).^2).*V;
    dVddx(V==0)  = 0;
    
    dVddy        = (-2*a+(2*a*y).^2).*V;
    dVddy(V==0)  = 0;
    
    dVdxdy       = 4*a*x.*y.*V;
    dVdxdy(V==0) = 0;
    
    
    %dVddx       = -2*tanh(z).*(1-(tanh(z).^2))/c^2;    
        
    %dVdth       = -dVdx.*r.*sin(th);
    %dVdr        = dVdx.*cos(th);
    
    %dVddth      = -dVdx.*r.*cos(th)+dVddx.*(r.*sin(th)).^2;
    %dVddr       = dVddx.*(cos(th)).^2;
    %dVdrdth     = -dVdx.*sin(th)-dVddx.*r.*sin(th).*cos(th);
    
    VDiff      = struct('dy1',dVdx,'dy2',dVdy,'ddy1',dVddx,'ddy2',dVddy,'dy1dy2',dVdxdy);%'Lap',dVddx
     
    VInt       = pi/(2*sqrt(a));
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

end