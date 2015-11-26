function [q,error] = HScorrection3N(dt,q0,v,D)

    
    r = rand(size(q0));
    
    ic = 0.1*ones(size(q0));
        
    q = fzero(@toMinimize,ic);
    
    error = cumPDF(q,dt,q0,v,D) - r;
        
    function F = toMinimize(q)
        F = cumPDF(q,dt,q0,v,D)-r;
    end

    function F = cumPDF(q,dt,q0,v,D)
        F = 1/2 * erfc( (q0 + v*dt) / sqrt(4*D*dt) ) ...
            - 1/2 * exp( v*q/D ) .* erfc( (q + q0 + v*dt) / sqrt(4*D*dt) );
        u = 1/2 * erfc( (q0 + v*dt) / sqrt(4*D*dt) );
        F = F./u;
    end

end