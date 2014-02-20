function w = MultMatrVecOp(A,v)
    [p,h1,h2]  = size(A);
    w        = zeros(p,p,2);

    w(:,:,1) = diag(A(:,1,1))*v(:,:,1) + diag(A(:,1,2))*v(:,:,2);
    w(:,:,2) = diag(A(:,2,1))*v(:,:,1) + diag(A(:,2,2))*v(:,:,2);
end