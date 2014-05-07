function K11 = JeffreyOnishi11Spherical(rf,rt,optsPhys)

    sigma   = optsPhys.sigmaS;
    sigmaH  = optsPhys.sigmaHS;
    lambda  = optsPhys.lambda;
    fx      = optsPhys.fx;
    fy      = optsPhys.fy;
    nMax    = optsPhys.nMax;


    r2m     = bsxfun(@minus, rf.^2 , rt.^2 );
    r2m2    = r2m.^2;

    % integration limits
    rm = bsxfun(@plus, rf , -rt );
    rp = bsxfun(@plus, rf ,  rt );

    lp=max(sigma^2,rp.^2);
    lm=max(sigma^2,rm.^2);
    
       % initialization
    K11=0;

    % use formulae from Jeffrey and Onishi evaluated at the end points
    for n=2:2:nMax
        K11=K11 + fx(n)*( X11(lp,n) - X11(lm,n) ) ...
            + fy(n)*( Y11(lp,n) - Y11(lm,n) );
    end

    %-------------------------------------------------------------------------
    % Spherically averaged integrals (see notes)
    %-------------------------------------------------------------------------
    
    function K=X11(u,k)
        if (k==2)
            K = (u - r2m2.*u.^(-1) + 2*r2m.*log(u) )/4;
        elseif (k==4)
            K = (-4*r2m.*u.^(-1) + 2*log(u) - r2m2.*u.^(-2) )/8;
        else
            K = -r2m2.*u.^(-k/2)/2/k -r2m.*u.^(-k/2+1)/(k-2) - u.^(-k/2+2)/(k-4)/2;
        end

        K = pi*(1+lambda)^(-k)*(sigmaH/2)^k*K.*bsxfun(@times,1./(rf.^3),rt);
    end


    function K=Y11(u,k)
        if (k==2)
            K=log(u);
        else
            K=-2*u.^(1-k/2)/(k-2); 
        end

        K=pi*(1+lambda)^(-k)*(sigmaH/2)^k*K.*bsxfun(@times,1./rf,rt);

        K = K - X11(u,k);
    end

    
end