function q = HSnew(dt,q0,v,D,N)
    
    q = q0 + v*dt + sqrt(2*D*dt) *randn(N,1);

    wallMask = (q<0);
    nCorrections = nnz(wallMask);
    
    qNew = zeros(nCorrections,1);
    errors = zeros(nCorrections,1);
    
    for iC = 1:nCorrections
%         [qTemp,error] = HScorrection(dt,q0,v,D);
%        [qTemp,error] = HScorrection2(dt,q0,v,D);
         [qTemp,error] = HScorrection3(dt,q0,v,D);
        qNew(iC) = qTemp;
        errors(iC) = error;
    end
        
    max(abs(errors))
    
    q(wallMask) = qNew;
    
end