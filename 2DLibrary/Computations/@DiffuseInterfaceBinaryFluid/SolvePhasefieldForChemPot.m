function mu = SolvePhasefieldForChemPot(this,uv,phi)

    m    = this.optsPhys.mobility;
    
    M    = this.IC.M;
    Diff = this.IC.Diff;
    Ind  = this.IC.Ind;
    
    rhs = (uv(1:M).*(Diff.Dy1*phi) + uv(M+1:end).*(Diff.Dy2*phi))/m;
    
    A   = Diff.Lap;
    
    %% BC1
    A(Ind.top|Ind.bottom,:) = 0;
	direction               = [zeros(M),eye(M)];
    A(Ind.top|Ind.bottom,:) = direction(Ind.top|Ind.bottom,:)*Diff.grad;    
    rhs(Ind.top|Ind.bottom) = 0;
    
    %% BC2
    EYM = eye(M);
    A(Ind.left|Ind.right,:) = EYM(Ind.left|Ind.right,:);
    rhs(Ind.left|Ind.right) = 0;           
    
    mu = A\rhs;

end