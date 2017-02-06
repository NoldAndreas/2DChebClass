function [D,DD] = SelfWallTermLS(x,y,optsPhys)

    sigmaHS = optsPhys.sigmaHS;
    
    yInv = y.^(-1);
    
    % no slip, keeping factor of D0 for main code
    % See Lauga & Squires, PoF, 17, 103102 (2005)
    
    Dx = (1-9/32*sigmaHS*yInv);
    Dy = (1-9/16*sigmaHS*yInv); % Note this becomes negative!!
    
    Dx_x = zeros(size(x));
    Dy_y = 9/16*sigmaHS*(yInv).^2;
    
    D  = [Dx;Dy];
    DD = [Dx_x;Dy_y]; 

end