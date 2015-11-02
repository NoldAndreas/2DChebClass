function theta = GHR_Inv(G,lambdaEta)

    options = optimoptions('fsolve','Display','off','TolFun',1e-8,'TolX',1e-8);    
    if(lambdaEta == 0)
        %theta = fsolve(@GHR_lambda0,ones(size(G)),options);
        theta = fsolve(@GHR,ones(size(G)),options);
    else
        theta = fsolve(@GHR,ones(size(G)),options);
    end
    
    %bisection(@GHR,0,pi,G);
    %fsolve(@GHR,pi/2*ones(size(G)),options);
    
    function z = GHR(t)
        z = GHR_lambdaEta(t,lambdaEta) - G;
    end

    function z = GHR_lambda0(t)
        z = GHR_lambdaEta0(t) - G;
    end

end