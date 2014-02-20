function K12 = RotnePrager12Spherical(r,s,optsPhys)

    alpha   = optsPhys.alpha;
    sigma   = optsPhys.sigma;
    sigmaH  = optsPhys.sigmaH;
    
    %------------------------------------------------------------------
    % set up frequently-used quantities
    %------------------------------------------------------------------

    r2m     = bsxfun(@minus, r.^2 , s.^2 );
    r2m2    = r2m.^2;
    r2p     = bsxfun(@plus, r.^2 , s.^2 );

    % integration limits
    rm = bsxfun(@plus, r , -s );
    rp = bsxfun(@plus, r ,  s );

    lp=max(sigma^2,rp.^2);
    lm=max(sigma^2,rm.^2);

    %------------------------------------------------------------------
    % evaluate integrand at endpoints
    %------------------------------------------------------------------


    K12 = getK12(lp)-getK12(lm);

    K12 = bsxfun(@times, pi./(r.^2), K12);

    function K12u=getK12(u)
       K12u = -3/16*sigmaH .* u.^(3/2) ...
              + ( 3/8*sigmaH*r2p + alpha/16*sigmaH^3 ) .* u.^(1/2) ...
              + ( -3/16*sigmaH*r2m2 - alpha/8*sigmaH^3*r2p ) .* u.^(-1/2) ...
              + alpha/16*sigmaH^3*r2m2 .* u.^(-3/2);
    end
        

    
end

