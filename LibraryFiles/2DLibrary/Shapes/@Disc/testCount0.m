 
    function c = testCount0(this,N) 
        s = struct('Ch',2*ones(N,1));                  
        c = 1;
        for i = 1:N
            c = c + s.Ch(i);
        end
    end   