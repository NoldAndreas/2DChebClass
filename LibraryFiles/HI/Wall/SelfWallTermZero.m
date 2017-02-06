function [D,DD] = SelfWallTermZero(x,y,optsPhys)
        
    Dx = ones(size(x));
    Dy = ones(size(x));
    
    Dx_x = zeros(size(x));
    Dy_y = zeros(size(x));
        
    D  = [Dx;Dy];
    DD = [Dx_x;Dy_y]; 

end