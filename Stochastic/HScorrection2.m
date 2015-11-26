function [q,error] = HScorrection2(dt,q0,v,D)

    dx = 1e-4;
    x = 0:dx:2;
    
    r = rand(1);

    F = toMinimize(x);
    
    [~,minPos] = min(abs(F));
    
    q = x(minPos);
    
    
%     ic = 0.1*ones(size(q0));
%         
%     fsolveOpts=optimset('Display','off');
%     
%     q = fsolve(@toMinimize,ic,fsolveOpts);
    
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