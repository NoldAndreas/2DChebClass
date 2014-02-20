function w = MultMatrVec(A,v)
    w = zeros(size(v));
    
    w(:,1) = A(:,1,1).*v(:,1) + A(:,1,2).*v(:,2);
    w(:,2) = A(:,2,1).*v(:,1) + A(:,2,2).*v(:,2);
end