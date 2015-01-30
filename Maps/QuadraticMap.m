function [z,dz] = QuadraticMap(x)        
    z  = (1/16)*x.*(5*x.^4-18*x.^2+29); 
    dz = (29/16)-(27/8)*x.^2+(25/16)*x.^4;
end