function [x] = InvExpMap(y,L,a)%dx,ddx,dddx,ddddx
 
    h    = L/(exp(1)-exp(1/4));
    z    = log(exp(1/4)+(y-a)/h);
    x    = 1-1./z;
%     h    = L/(exp(1)-exp(1/2));    
%     z    = log(exp(1/2)+(y-a)/h);
%     x    = 1-1./z;        
        
end