function [q,error] = HSCollisionSampling(dt,q0,v,D)

    
    r = rand(1);
    
    ic = 0.1;
        
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