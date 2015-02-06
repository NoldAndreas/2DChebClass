function theta = f_0_DussanRameGaroff(r,thOut,alpha,R_T)
    % r and R_T are given in multiples of the capillary length
    beta  = alpha - thOut;
    theta = thOut + ...
            r*(2*sin(beta/2)+ 1/R_T*( 2*((cos(beta/2))^3-1)/(3*sin(beta/2)) +sin(beta) )  ) ...
            - r.^2*sin(beta)/2;    
end
