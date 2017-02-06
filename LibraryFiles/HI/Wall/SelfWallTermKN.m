function [D,DD] = SelfWallTermKN(x,y,optsPhys)

    sigmaHS = optsPhys.sigmaHS;
    
    a = sigmaHS/2;
    
    yInv = y.^(-1);
    
    % no slip, keeping factor of D0 for main code
    % See vonHansen, Hinczewski, Netz, JCP 134, 235102 (2011)
    
    Dx = (1 - 9/16*a*yInv + 1/8*(a*yInv).^3);
    Dy = (1 - 9/8*a*yInv + 1/2*(a*yInv).^3); 
    
    Dx_x = zeros(size(x));
    Dy_y = 9/8*a*(yInv).^2 - 3/2*a^3*yInv.^4;
    
    D  = [Dx;Dy];
    DD = [Dx_x;Dy_y]; 


end