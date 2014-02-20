function [z,dz,dx,ddx,dddx] = MapInf2(x,L)
%map
       z    =  L*x./(1-x.^2);
       dz   =  L*(1+x.^2)./((1-x.^2).^2);
       if (nargout == 2), return; end
       h    = sqrt(L^2 + 4*z.^2);
       dx   = L./(2*z.^2).*(1 - L./h);
       ddx  = L./(z.^3).*(-1 + L.*(L^2 + 6*z.^2)./(h.^3));
       dddx = 3*(L-h)./(z.^4) + 6./(h.*z.^2) + 96*(z.^2)./(h.^5);

end