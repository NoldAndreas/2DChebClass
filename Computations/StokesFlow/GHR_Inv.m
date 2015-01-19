function theta = GHR_Inv(G,lambda)

    theta = fsolve(@GHR,pi/2*ones(size(G)));
    
    function z = GHR(t)
        z = GHR_lambda(t,lambda) - G;
    end

end