function D = DiffusionCoefficientWall(x,y,optsPhys)

    a = optsPhys.sigmaHS/2;
    
    wallPos = optsPhys.wallPos;
    
    h = y - wallPos;
    hInv = h.^(-1);

    % no slip, keeping factor of D0 for main code
    
    Dy = (1-9/8*a*hInv);
    Dx = (1-9/16*a*hInv);
    
    D = [Dx;Dy];
end