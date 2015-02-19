function [z,dzdr_r,alpha] = ConstShortRange(r,parameter) 
    z      = ones(size(r));
    dzdr_r = zeros(size(r));
    
    alpha  = -2/3*pi*(parameter.lambda^3 - parameter.LJsigma^3);
    
end