function [y,Diff] = M1Tref(x,d,ep)

    if(length(ep) == 1)
        ep = ep*ones(size(x));
    end
    mark = (ep == inf);
        
    c1  = (asinh((1-d)./ep)+asinh((1+d)./ep))/2;
    c2  = asinh((1-d)./ep);

    y          = d + ep.*sinh(c1.*(x-1)+c2);        
    y(mark)    = x(mark);
    y(x==1)    = 1;
    y(x==-1)   = -1;

    Diff.dydx      = ep.*c1.*cosh(c1.*(x-1)+c2);
    Diff.dydx(mark) = 1;

    Diff.dyddx     = ep.*(c1.^2).*sinh(c1.*(x-1)+c2);
    Diff.dyddx(mark) = 0;

    dh1        = (1./(sqrt(ep.^2+(1+d).^2)) - 1./(sqrt(ep.^2+(1-d).^2)))/2;
    dh2        = 1./(sqrt(ep.^2+(1-d).^2));

    de1               = -((1-d)./((ep.^2+(1-d).^2).^(3/2)) + (1+d)./((ep.^2+(1+d).^2).^(3/2)))/2;
    de2               = (1-d)./((ep.^2+(1-d).^2).^(3/2));

    Diff.dyd_d        = 1 + ep.*cosh(c1.*(x-1)+c2).*((x-1).*dh1-dh2);
    Diff.dyd_d(mark)  = 0;
    Diff.dydd_d       = ep.*sinh(c1.*(x-1)+c2).*((x-1).*dh1-dh2).^2 + ...
                        ep.*cosh(c1.*(x-1)+c2).*((x-1).*de1-de2);             
    Diff.dydd_d(mark) = 0;                     

    Diff.dydxd_d      = ep.*c1.*sinh(c1.*(x-1)+c2).*(dh1.*(x-1)-dh2) +...
                        ep.*dh1.*cosh(c1.*(x-1)+c2);
    Diff.dydxd_d(mark) = 0;

    de1               = -((1-d)./(ep.*sqrt(ep.^2+(1-d).^2)) + (1+d)./(ep.*sqrt(ep.^2+(1+d).^2)))/2;
    de2               = (1-d)./(ep.*sqrt(ep.^2+(1-d).^2));

    z                 = c1.*(x-1)+c2;
    Diff.dyde         = sinh(z) + ep.*cosh(z).*(de1.*(x-1)-de2);
    Diff.dyde(mark)   = 0;
    
    Diff.dydedx       = cosh(z).*c1+ep.*sinh(z).*c1.*((x-1).*de1-de2) + ep.*cosh(z).*de1;
    Diff.dydedx(mark) = 0;
    
    
    df1               = (1-d)./(ep.^2.*(ep.^2+(1-d).^2).^(1/2)) - 1/2*(1-d).^3./(ep.^2.*(ep.^2+(1-d).^2).^(3/2)) + ...
                        +(1+d)./(ep.^2.*(ep.^2+(1+d).^2).^(1/2)) - 1/2*(1+d).^3./(ep.^2.*(ep.^2+(1+d).^2).^(3/2));

    df2               = 2*(1-d)./(ep.^2.*(ep.^2+(1-d).^2).^(1/2)) - (1-d).^3./(ep.^2.*(ep.^2+(1-d).^2).^(3/2));
    
    Diff.dydde        = 2*cosh(z).*((x-1).*de1-de2)+...
                        +ep.*sinh(z).*((x-1).*de1-de2).^2+...
                        +ep.*cosh(z).*((x-1).*df1+df2);
    Diff.dydde(mark)  = 0;
    
    %***************************************
    %diff(dy,d,ep):
    Z     = c1.*(x-1)+c2;
    DM    = 1./sqrt(1+(1-d).^2./ep.^2);
    DP	  = 1./sqrt(1+(1+d).^2./ep.^2);
    dh1   = (DP-DM)./(2*ep);
    dh2   = DM./ep;
    ddh1  = -(DM.*(1-d)+DP.*(1+d))./(2*ep); 
    ddh2  = DM.*(1-d)./ep; 
    dddh1 = (DM-DM.^3.*(1-d).^2./ep.^2-DP+DP.^3.*(1+d).^2./ep.^2)./(2*ep);
    dddh2 = -DM./ep+DM.^3.*(1-d).^2./ep.^3;
    H1    = dh1.*(x-1)-dh2; 
    H2    = ddh1.*(x-1)-ddh2; 
    H3    = dddh1.*(x-1)-dddh2;
    
    Diff.dyd_ddep = cosh(Z).*(H1+H3)+sinh(Z).*H2.*H1;
    %***************************************
        
end