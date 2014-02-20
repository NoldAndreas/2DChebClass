function Diff = PhysicalDerivatives_1D(PhysSpace,x,CompDiff)    

    Diff = struct();

    [h_1,h_2,dx1,ddx1,dddx1,ddddx1] = PhysSpace(x); 
    Diff.Dy   = sparse(diag(dx1)*CompDiff.Dx);  
    
    if(isfield(CompDiff,'DDx'))    
        Diff.DDy      = sparse(diag(ddx1)*CompDiff.Dx + (diag(dx1.^2))*CompDiff.DDx);    
    end
    if(isfield(CompDiff,'DDDx'))
        Diff.DDDy     = diag(dx1.^3)*CompDiff.DDDx + 3*diag(dx1.*ddx1)*CompDiff.DDx + diag(dddx1)*CompDiff.Dx;       
    end
    if(isfield(CompDiff,'DDDDx'))        
        Diff.DDDDy    = diag(dx1.^4)*CompDiff.DDDDx + ...
                        6*diag(dx1.^2.*ddx1)*CompDiff.DDDx + ...
                        3*diag(ddx1.^2)*CompDiff.DDx + ...
                        4*diag(dx1.*dddx1)*CompDiff.DDx + ...
                        diag(ddddx1)*CompDiff.Dx;
    end
end
