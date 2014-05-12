function [z,dzdr_r] = Gaussian(r,optsPhys)
    %z      - the 2-particle interaction, distance r
    %dzdr_r - 1/r * dz/dr
    epsilon = optsPhys.epsilon;
    alpha   = optsPhys.alpha;
    
	z = epsilon.*exp(-(r.^2)./(alpha.^2));
    z(r == inf) = 0;
        
    dzdr_r = -2*epsilon./(alpha.^2).*exp(-(r.^2)./(alpha.^2));  
end
