function z = EricV2(x,y,optsPhys)

    V0      = optsPhys.V0;
    sigma   = optsPhys.sigma;
    
	% z = epsilon*exp(-(r.^2)/(alpha^2)) but need to do spherical averaging
    % we do epsilon \int y^2 exp(-1/alpha^2 (x^2 + y^2 - 2 x y cos(t))
    % sin(t) dt dp:
    
    z = - pi * V0 * sigma^2 .* y./x .* ...
          ( exp(-(x+y).^2/sigma^2) - exp(-(x-y).^2/sigma^2) );
    
    z(x == inf | y == inf) = 0;
    
end