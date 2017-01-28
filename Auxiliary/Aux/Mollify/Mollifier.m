function z = Mollifier(r)        
        z       = zeros(size(r));
        z(r<1)  = exp(-1./(1-(r(r<1)).^2));
        z(r>=1) = 0;
end