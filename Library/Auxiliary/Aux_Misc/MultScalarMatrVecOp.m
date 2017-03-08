function w = MultScalarMatrVecOp(A,vOp)    
    w        = zeros(size(vOp));

    w(:,:,1) = A(1,1)*vOp(:,:,1) + A(1,2)*vOp(:,:,2);
    w(:,:,2) = A(2,1)*vOp(:,:,1) + A(2,2)*vOp(:,:,2);
end