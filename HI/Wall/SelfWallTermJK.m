function [D,DD] = SelfWallTermJK(x,y,optsPhys)

    sigmaHS = optsPhys.sigmaHS;
    
    yInv = y.^(-1);
    
    % no slip, keeping factor of D0 for main code
    % See Jones & Kutteh, evaluate additional term for self interaction,
    % i.e. with r'=r and ignoring first Oseen term
    
    Dx = (1-3/32*sigmaHS*yInv);
    Dy = (1-3/8*sigmaHS*yInv);
    
    Dx_x = zeros(size(x));
    Dy_y = 3/8*sigmaHS*(yInv).^2;
    
    D  = [Dx;Dy];
    DD = [Dx_x;Dy_y]; 

end