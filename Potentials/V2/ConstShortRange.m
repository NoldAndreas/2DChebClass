function [z,dzdr_r,alpha] = ConstShortRange(r,V2) 
    epsilon = V2.epsilon;

    lambda  = 1.5;
    LJsigma = 1;       

    z      = ones(size(r));
    dzdr_r = zeros(size(r));        
    alpha  = -2/3*pi*(lambda^3 - LJsigma^3);    
    
    c      = epsilon*(-16/9*pi)/alpha;    
    z      = c*z;
    dzdr_r = c*dzdr_r;
    alpha  = c*alpha;
    
end