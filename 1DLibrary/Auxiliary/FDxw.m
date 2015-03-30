function [x,w,D] = FDxw(N)
    
    x=(-1:2/(N-1):1)'; 

    h=2/(N-1);

    w=[h/2 h*ones(1,N-2) h/2];
    
    D = sparse(2:N-1,[(3:N)],1/(2*h),N,N) ...
            - sparse(2:N-1,[1:(N-2)],1/(2*h),N,N);
        
	D(1,1) = -1/h; D(1,2) = 1/h;
    D(N,N) = 1/h;  D(N,N-1) = -1/h;
        
end