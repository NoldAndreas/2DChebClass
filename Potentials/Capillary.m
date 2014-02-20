function V = Capillary(h,epsilon_w)
%V = Capillary(h,epsilon_w)
% For capillary of width 2*h
    alpha = - epsilon_w*pi^2/2;    
    V     = 2*alpha - 2*AttractiveWall(h,epsilon_w,1);
end