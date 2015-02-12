function x = InvLogMap(y,L,a)%dx,ddx,dddx,ddddx
    
    h = L/log(1/2);
    x = 1 - 2*exp((y-a)/h);
        
end