function u = BH_1D_I(y,r)
    u = 2*pi*( 1./(y.^2 + r.^2).^2 - 1./(5*(y.^2 + r.^2).^5) );
end