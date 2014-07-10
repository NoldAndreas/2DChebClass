function [D,DD] = DiffusionCoefficientWall(x,y,optsPhys)

    a = optsPhys.sigmaHS/2;
    
    wallPos = optsPhys.wallPos;
    
    h = y - wallPos;
    hInv = h.^(-1);

    % no slip, keeping factor of D0 for main code
    % Lauga & Squires, PoF 17, 103102 (2005)
    
    Dx = (1-9/16*a*hInv);
    Dy = (1-9/8*a*hInv);
    
    Dx_x = zeros(size(x));
    Dy_y = 9/8*a*(hInv).^2;
    
    D  = [Dx;Dy];
    DD = [Dx_x;Dy_y]; 

end