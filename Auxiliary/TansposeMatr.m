function C = TansposeMatr(A)
    [p,h1,h2] = size(A);
    C = zeros(p,2,2);
    
    C(:,1,1) = A(:,1,1);
    C(:,1,2) = A(:,2,1);
    C(:,2,1) = A(:,1,2);
    C(:,2,2) = A(:,2,2);
end