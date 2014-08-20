function DiffOut = PhysicalDerivatives(a_shape,Pts,Sel,CompDiff1,CompDiff2,order)%Dx1,Dx2,DDx1,DDx2,DDDx1,DDDx2)

    n1       = length(Pts.x1);
    n2       = length(Pts.x2);
    
    Diff1 = PhysicalDerivatives_1D(@a_shape.PhysSpace1,Pts.x1,CompDiff1);
    Diff2 = PhysicalDerivatives_1D(@a_shape.PhysSpace2,Pts.x2,CompDiff2);
        
    Diff = struct();
    
    % conversion into physical space:     
    %dy1       = sparse(diag(dx1)*Dx1);   
    %dy2       = sparse(diag(dx2)*Dx2); 
    Diff.Dy2  = sparse(kron(eye(n1),Diff2.Dy));   
    Diff.Dy1  = sparse(kron(Diff1.Dy,eye(n2))); 
    Diff.grad = sparse([Diff.Dy1;Diff.Dy2]);            
    Diff.div  = sparse([Diff.Dy1 Diff.Dy2]);
    
    if (order == 1)        
        CopyFields();
        return;
    end
    
    %ddy1    = sparse(diag(ddx1)*Dx1 + (diag(dx1.^2))*DDx1);    
    %ddy2    = sparse(diag(ddx2)*Dx2 + (diag(dx2.^2))*DDx2);        
    
    Diff.DDy1     = sparse(kron(Diff1.DDy,eye(n2)));
    Diff.DDy2     = sparse(kron(eye(n1),Diff2.DDy));            
    Diff.Dy1Dy2   = kron(Diff1.Dy,Diff2.Dy);    
    Diff.Lap      = sparse(Diff.DDy1 + Diff.DDy2);    
    Diff.gradDiv  = [Diff.DDy1 Diff.Dy1Dy2 ; Diff.Dy1Dy2 Diff.DDy2];
    Diff.LapVec   = [Diff.Lap zeros(n1*n2) ; zeros(n1*n2) Diff.Lap];
       
    if (order == 2)        
        CopyFields(); 
        return;
    end
    
    %dddy1        = diag(dx1.^3)*DDDx1 + 3*diag(dx1.*ddx1)*DDx1 + diag(dddx1)*Dx1;       
    %dddy2        = diag(dx2.^3)*DDDx2 + 3*diag(dx2.*ddx2)*DDx2 + diag(dddx2)*Dx2;       
    
    Diff.Dy1DDy2 = kron(Diff1.Dy,Diff2.DDy);
    Diff.DDy1Dy2 = kron(Diff1.DDy,Diff2.Dy);    
    Diff.DDDy1   = kron(Diff1.DDDy,eye(n2));
    Diff.DDDy2   = kron(eye(n1),Diff2.DDDy);
    Diff.gradLap = [Diff.DDDy1 + Diff.Dy1DDy2 ; Diff.DDy1Dy2 + Diff.DDDy2];
    Diff.LapDiv  = [Diff.DDDy1 + Diff.Dy1DDy2 , Diff.DDy1Dy2 + Diff.DDDy2];
    
    if (order == 3)        
        CopyFields(); 
        return;
    end
    
    Diff.DDDDy1   = kron(Diff1.DDDDy,eye(n2));
    Diff.DDDDy2   = kron(eye(n1),Diff2.DDDDy);
    Diff.DDy1DDy2 = kron(Diff1.DDy,Diff2.DDy);
    Diff.Dy1DDDy2 = kron(Diff1.Dy,Diff2.DDDy);
    Diff.DDDy1Dy2 = kron(Diff1.DDDy,Diff2.Dy);
    Diff.Lap2     = Diff.DDDDy1 + 2*Diff.DDy1DDy2 + Diff.DDDDy2;
    
    CopyFields();
    
    function CopyFields()        
        DiffOut = struct();
        for i = 1:length(Sel)
            DiffOut.(Sel{i}) = Diff.(Sel{i});
        end                               
    end

end