function theta = GHR_Inv(G,lambdaEta)

    options = optimoptions('fsolve','Display','off'); 
    theta = fsolve(@GHR,pi/2*ones(size(G)),options);
    
    function z = GHR(t)
        z = GHR_lambdaEta(t,lambdaEta) - G;
    end

end