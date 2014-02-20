function V = Rectangle(x,y,epsilon)
%V = Rectangle(x,y,epsilon)
% V = int( phi_attr(xs,ys), xs = -x..x,ys = -y..y)
% where  phi_attr = Phi2DLongRange.m 

    V = -2*pi*epsilon*( atan(x.*y./(sqrt(1+x.^2+y.^2))) +...
                + x.*y.*(2+x.^2+y.^2)./((1+x.^2).*(1+y.^2).*sqrt(1+x.^2+y.^2)));
    
end