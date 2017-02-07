function [y,dy] = AttractiveWall(z,epsilon_w,sigma_w)
    z = z/sigma_w;
    y  = epsilon_w*pi*(atan(z)-pi/2+z./(1+z.^2));
    dy = epsilon_w*pi*2./((1+z.^2).^2);
end