function theta = GHR_Inv(G,lambdaEta)

    options = optimoptions('fsolve','Display','off','TolFun',1e-8,'TolX',1e-8);    
    theta = fsolve(@GHR,ones(size(G)),options);
    
    %bisection(@GHR,0,pi,G);
    %fsolve(@GHR,pi/2*ones(size(G)),options);
    
    function z = GHR(t)
        z = GHR_lambdaEta(t,lambdaEta) - G;
    end

end