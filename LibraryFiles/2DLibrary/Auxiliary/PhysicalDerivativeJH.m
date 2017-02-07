function Diff = PhysicalDerivativeJH(Pts,Diff1,Diff2,J,dH1,dH2)

    markInf  = (Pts.y1_kv == inf) | (Pts.y1_kv == - inf) ...
             | (Pts.y2_kv == inf) | (Pts.y2_kv == - inf);
         
    mark1Inf  = (Pts.y1_kv == inf) | (Pts.y1_kv == - inf);
    mark2Inf  = (Pts.y2_kv == inf) | (Pts.y2_kv == - inf);
    
    N1 = length(Diff1.Dx);    N2 = length(Diff2.Dx);

    Dx1             = sparse(kron(Diff1.Dx,eye(N2))); 
    Dx2             = sparse(kron(eye(N1),Diff2.Dx)); 
    
    DDx1            = sparse(kron(Diff1.DDx,eye(N2))); 
    DDx2            = sparse(kron(eye(N1),Diff2.DDx)); 
    Dx1Dx2          = sparse(kron(Diff1.Dx,Diff2.Dx)); 
        
    gradx           = zeros(N1*N2,N1*N2,2);
    gradx(:,:,1)    = Dx1;
    gradx(:,:,2)    = Dx2;
    
    JT              = TansposeMatr(J);
    JT_I            = InvMatr(JT);
    J_I             = InvMatr(J);
    
    dH1(mark1Inf,:,:)  = 0;
    dH2(mark2Inf,:,:)  = 0;
    JT_I(markInf,:,:) = 0;
    JT_I(markInf,1,1) = 1./JT(markInf,1,1);
    JT_I(markInf,2,2) = 1./JT(markInf,2,2);
    
    J_I(markInf,:,:)  = 0;
    J_I(markInf,1,1) = 1./J(markInf,1,1);
    J_I(markInf,2,2) = 1./J(markInf,2,2);

    
    gradh           = MultMatrVecOp(JT_I,gradx);
    Diff.grad       = zeros(2*N1*N2,N1*N2);
    Diff.grad(1:N1*N2,:)     = gradh(:,:,1);
    Diff.grad(1+N1*N2:end,:) = gradh(:,:,2);
    
    Diff.Dy1        = gradh(:,:,1);
    Diff.Dy2        = gradh(:,:,2);
    
    Hx = zeros(N1*N2,N1*N2,2,2);
    Hx(:,:,1,1) = DDx1;
    Hx(:,:,1,2) = Dx1Dx2;
    Hx(:,:,2,1) = Dx1Dx2;
    Hx(:,:,2,2) = DDx2;
        
    %gradh(markInf,markInf,1) = 0;
    
    A          = zeros(N1*N2,N1*N2,2,2);
    A(:,:,1,1) = diag(dH1(:,1,1))*gradh(:,:,1) + diag(dH2(:,1,1))*gradh(:,:,2);
    A(:,:,1,2) = diag(dH1(:,1,2))*gradh(:,:,1) + diag(dH2(:,1,2))*gradh(:,:,2);
    A(:,:,2,1) = diag(dH1(:,2,1))*gradh(:,:,1) + diag(dH2(:,2,1))*gradh(:,:,2);
    A(:,:,2,2) = diag(dH1(:,2,2))*gradh(:,:,1) + diag(dH2(:,2,2))*gradh(:,:,2);    
    
    Hy = Mult2MatrMatrOp(JT_I,Hx - A,J_I);
    Diff.DDy1   = Hy(:,:,1,1);
    Diff.Dy1Dy2 = Hy(:,:,1,2);
    Diff.DDy2   = Hy(:,:,2,2);
            
    %Diff.DDy1(markInf,:,:) = 0;
    
    Diff.Lap    = Hy(:,:,1,1) + Hy(:,:,2,2);
    Diff.div    = [Diff.Dy1 , Diff.Dy2];
    Diff.LapVec = [Diff.Lap , zeros(N1*N2); zeros(N1*N2) , Diff.Lap];
    Diff.gradDivVec = [Diff.DDy1 , Diff.Dy1Dy2 ; Diff.Dy1Dy2 , Diff.DDy2];

end