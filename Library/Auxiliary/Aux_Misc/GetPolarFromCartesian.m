function [ur utheta] = GetPolarFromCartesian(ux,uy,theta)

% y
% I
% I
% I
% I------------> x
% I
% I
% I

    ur     = ux.*cos(theta)+uy.*sin(theta);
    utheta = -ux.*sin(theta)+uy.*cos(theta);

end