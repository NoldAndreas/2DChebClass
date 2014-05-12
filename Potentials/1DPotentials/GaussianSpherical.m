function z = GaussianSpherical(x,y,optsPhys)

    epsilon = optsPhys.epsilon;
    alpha   = optsPhys.alpha;
    
	% z = epsilon*exp(-(r.^2)/(alpha^2)) but need to to spherical averaging
    % we do epsilon \int y^2 exp(-1/alpha^2 (x^2 + y^2 - 2 x y cos(t))
    % sin(t) dt dp:
    
    z = - pi * epsilon * alpha^2 .* y./x .* ...
          ( exp(-(x+y).^2/alpha^2) - exp(-(x-y).^2/alpha^2) );
    
    z(x == inf | y == inf) = 0;
    
end