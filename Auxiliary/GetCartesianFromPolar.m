function [u,v] = GetCartesianFromPolar(ur,utheta,phi)

    u = ur.*cos(phi)-utheta.*sin(phi);
    v = ur.*sin(phi)+utheta.*cos(phi);

end