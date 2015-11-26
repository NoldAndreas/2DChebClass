function q = HSrej(dt,q0,v,D,N)
    
    q = q0 + v*dt + sqrt(2*D*dt) *randn(N,1);

    nnz(q<0)/N;
    
    q(q<0) = q0;
    

    
end