function w = MultScalarMatrOp(A,vOp)    
    p1       = size(vOp,1);
    p2       = size(vOp,3);
    w        = zeros(p1,p1,p2,p2);

    w(:,:,1,1) = A(1,1)*vOp(:,:,1,1) + A(1,2)*vOp(:,:,2,1);    
    w(:,:,1,2) = A(1,1)*vOp(:,:,1,2) + A(1,2)*vOp(:,:,2,2);    
    w(:,:,2,1) = A(2,1)*vOp(:,:,1,1) + A(2,2)*vOp(:,:,2,1);    
    w(:,:,2,2) = A(2,1)*vOp(:,:,1,2) + A(2,2)*vOp(:,:,2,2);
end