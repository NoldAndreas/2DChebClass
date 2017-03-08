function K12 = JeffreyOnishi12Spherical(rf,rt,optsPhys)

    sigma   = optsPhys.sigmaS;
    sigmaH  = optsPhys.sigmaHS;
    lambda  = optsPhys.lambda;
    fx      = optsPhys.fx;
    fy      = optsPhys.fy;
    nMax    = optsPhys.nMax;

    r2m     = bsxfun(@minus, rf.^2 , rt.^2 );
    r2m2    = r2m.^2;
    r2p     = bsxfun(@plus, rf.^2 , rt.^2 );


    % integration limits
    rm = bsxfun(@plus, rf , -rt );
    rp = bsxfun(@plus, rf ,  rt );

    lp=max(sigma^2,rp.^2);
    lm=max(sigma^2,rm.^2);

    
    % initialization
    K12=0;

    % use formulae from Jeffrey and Onishi evaluated at the end points
    for n=1:2:nMax
        K12=K12 - fx(n)*( X12(lp,n) - X12(lm,n) ) ...
            - fy(n)*( Y12(lp,n) - Y12(lm,n) );
    end

    %-------------------------------------------------------------------------
    % Spherically averaged integrals (see notes)
    %-------------------------------------------------------------------------

    function K=X12(u,k)
        if (k==4)
            K = -r2m2.*u.^(-2)/2 - log(u);
        else
            K = -2*r2m2.*u.^(-k/2)/k + 2*u.^(-k/2+2)/(k-4);
        end

        K = pi/4*(1+lambda)^(-k)*(sigmaH/2)^k*K.*bsxfun(@times,1./(rf.^2),ones(size(rt)));
    end


    function K=Y12(u,k)
        if (k==2)
            K = 2*r2p.*log(u) - 2*u;
        elseif (k==4)
            K = -2*r2p.*u.^(-1) - 2*log(u);
        else
            K = -4*r2p.*u.^(1-k/2)/(k-2) + 4*u.^(2-k/2)/(k-4);
        end

        K=pi/4*(1+lambda)^(-k)*(sigmaH/2)^k*K.*bsxfun(@times,1./(rf.^2),ones(size(rt)));

        K=K - X12(u,k);
    end

    
end