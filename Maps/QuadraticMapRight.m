function [z,dz] = QuadraticMapRight(x)
        % z = L*((xf+1)/2)^2   
    z  = (1 - (((1-x)/2).^2));%- (((1+xf)/2).^2));
    dz = (1-x)/2;
   % dz = (29/16)-(27/8)*x.^2+(25/16)*x.^4;
end