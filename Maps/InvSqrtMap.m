function x = InvSqrtMap(z,L05,LD)
%[-LD,LD] -> [-1,1]       and   
%hS12     -> sqrt(0.5)
    h2   = L05^2/(  1 - 2*L05^2/LD^2);
    C    = sqrt(1+h2/LD^2);
    x    = C*z.*(h2+z.^2).^(-1/2);
    
    x(z==-Inf) = -C;
    x(z==Inf)  = C;

end