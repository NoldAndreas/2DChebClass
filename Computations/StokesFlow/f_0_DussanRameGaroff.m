function theta = f_0_DussanRameGaroff(r,thOut,alpha,R_T)
    % r and R_T are given in multiples of the capillary length
    beta  = alpha - thOut;
    theta = thOut + 2*r*sin(beta/2) - r.^2*sin(beta)/2 ...
            + 1/R_T*( r.*( 2/3*((cos(beta/2))^3-1)/sin(beta/2) +sin(beta) ) ...
                     -r.^2*(sin(beta/2)*cos(beta)));    
end
