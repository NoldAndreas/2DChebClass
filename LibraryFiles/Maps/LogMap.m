function [y,dy] = LogMap(x,L,a)%dx,ddx,dddx,ddddx
% z = a + ( L*exp(1/(1-x)-1) - exp(-1/2) )
    
    h = L/log(1/2);
    y = a + h*log((1-x)/2);
    dy = h./(x-1);    
        
end