function [z,dzdr_r,a] = Gaussian(r,optsPhys)
    %z      - the 2-particle interaction at distance r

    %%
    % 
    % $$z = \varepsilon e^{-r^2/\alpha^2}$$
    % 
    % $$dzdr_r = \frac{1}{r} \frac{dz}{dr}$$
    %
    % $$ a  = \pi \int_0^\infty  r z(r) dr = \frac{\alpha^2}{2} \pi \varepsilon$$
    %
    
    
    epsilon = optsPhys.epsilon;
    alpha   = optsPhys.alpha;
    
	z = epsilon.*exp(-(r.^2)./(alpha.^2));
    z(r == inf) = 0;
        
    dzdr_r = -2*epsilon./(alpha.^2).*exp(-(r.^2)./(alpha.^2));  
    
    a = pi*epsilon*alpha^2/2;
end
