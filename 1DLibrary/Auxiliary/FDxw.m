function [x,w] = FDxw(N)
    
    x=(-1:2/(N-1):1)'; 

    h=2/(N-1);

    w=[h/2 h*ones(1,N-2) h/2];
end