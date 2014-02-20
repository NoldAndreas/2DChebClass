function [u,v] = GetCartesianFromPolar(ur,utheta,theta)

% y
% I
% I
% I
% I------------> x
% I
% I
% I

    u = ur.*cos(theta)-utheta.*sin(theta);
    v = ur.*sin(theta)+utheta.*cos(theta);

end