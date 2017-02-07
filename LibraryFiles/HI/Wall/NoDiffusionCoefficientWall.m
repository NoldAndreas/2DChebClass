function [D,DD] = NoDiffusionCoefficientWall(x,y,optsPhys)

    % no wall effects - just for testing accuracy compared to standard
    % calculation

    Dx = ones(size(x));
    Dy = ones(size(y));

    Dx_x = zeros(size(x));
    Dy_y = zeros(size(y));
    
    D  = [Dx;Dy];
    DD = [Dx_x;Dy_y]; 
end