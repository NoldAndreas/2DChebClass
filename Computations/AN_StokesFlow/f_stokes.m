function [f,dfdt] = f_stokes(t,lambda)

    S = sin(t);
    C = cos(t);
   
    A        = 2*sin(t).*(lambda^2*(t.^2-(sin(t)).^2)+2*lambda*(t.*(pi-t)+(sin(t)).^2)+(pi-t).^2-(sin(t)).^2);
    dA       = 2*cos(t).*(lambda^2*(t.^2-(sin(t)).^2)+2*lambda*(t.*(pi-t)+(sin(t)).^2)+(pi-t).^2-(sin(t)).^2) + ...
               2*sin(t).*(lambda^2*(2*t-2*S.*C) + 2*lambda*(pi-2*t+2*S.*C)-2*(pi-t)-2*S.*C);
    
    B        = lambda*(t.^2-(sin(t)).^2).*(pi-t+sin(t).*cos(t))+((pi-t).^2-(sin(t)).^2).*(t-sin(t).*cos(t));
    dB       = lambda*(2*t-2*S.*C).*(pi-t+sin(t).*cos(t))+ ...
               lambda*(t.^2-(sin(t)).^2).*(-1+C.^2-S.^2) +...
               (-2*(pi-t)-2*S.*C).*(t-sin(t).*cos(t))    +...
               ((pi-t).^2-(sin(t)).^2).*(1-C.^2+S.^2);

    f        = A./B;              
    dfdt     = dA./B - A.*dB./(B.^2);
    
    
    f(t==0)  = inf;
    f(t==pi) = inf;
end