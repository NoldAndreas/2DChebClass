function det = DetMatr(A)        
    det     = A(:,1,1).*A(:,2,2)-A(:,1,2).*A(:,2,1);
end