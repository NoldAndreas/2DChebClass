function [wall,dwall] = BarkerHendersonWall(y2,epsilon_w)
    
    y2             = y2 + 1;
    wall           = 4*pi*(1./(45*y2.^9) - 1./(6*y2.^3));
    wall(y2==inf)  = 0;    %wall(yw==-inf) = -pi^2;
    
    dwall          = -4*pi.*(1./(5*y2.^10)-1./(2*y2.^4));
    dwall(y2==inf) = 0;
    
    wall           = epsilon_w * wall;
    dwall          = epsilon_w * dwall;

end