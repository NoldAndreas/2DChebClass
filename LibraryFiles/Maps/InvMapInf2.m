function [x] = InvMapInf2(z,L)
        x            = (sqrt(L^2+4*z.^2)-L)./(2*z);
        x(z == inf)  = 1;
        x(z == -inf) = -1;
        x(z == 0)    = 0;
end