function [q,error] = newHSTest(dt,q0,v,D)

%     u = upperBound(dt,q0,v,D);
%     l = 0;
%    
%     r = l + (u-l)*rand(1);

    r = rand(1);

    ic = 0.1;
        
    fsolveOpts=optimset('Display','off');
    
    q = fsolve(@toMinimize,ic,fsolveOpts);
    
    error = cumPDF(q,dt,q0,v,D) - r;
        
    function F = toMinimize(q)
        F = cumPDF(q,dt,q0,v,D)-r;
    end

    function F = cumPDF(q,dt,q0,v,D)
        F = 1/2 * erfc( (q0 + v*dt) / sqrt(4*D*dt) ) ...
            - 1/2 * exp( v*q/D ) .* erfc( (q + q0 + v*dt) / sqrt(4*D*dt) );
        u = 1/2 * erfc( (q0 + v*dt) / sqrt(4*D*dt) );
        F = F/u;
    end

end