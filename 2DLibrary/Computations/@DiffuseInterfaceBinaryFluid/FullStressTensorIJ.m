function [A,b] = FullStressTensorIJ(this,phi,i,j)
    % get matrices for
    % T = Cak*(grad(u) + grad(u)^T) +...
    %         + (-p*Cak + W(phi) + Cn/2*|grad(phi)|^2 )*I - Cn*(grad(phi) X grad(phi))
    %   = A*[p;uv] + b

    Cn    = this.optsPhys.Cn;
    Cak   = this.optsPhys.Cak;    
    Diff  = this.IC.Diff;
    M     = this.IC.M;

    [~,W]   = DoublewellPotential(phi,Cn);
    bDiag   = W + Cn/2*((Diff.Dy1*phi).^2 + (Diff.Dy2*phi).^2); %CahnHilliardFreeEnergy(rho,Cn,Diff);    
    if(i == 1 && j == 1)
        Ap      = -Cak*eye(M);
        Auv     = 2*Cak*[Diff.Dy1 , zeros(M)];
        b       = bDiag - Cn*(Diff.Dy1*phi).^2;
    elseif((i==2 && j == 1)  || (i==1 && j == 2))
        Ap      = zeros(M);
        Auv     = Cak*[Diff.Dy2 , Diff.Dy1];
        b       = - Cn*(Diff.Dy1*phi).*(Diff.Dy2*phi);
    elseif(i==2 && j == 2)
        Ap      = -eye(M);
        Auv     = 2*Cak*[zeros(M) , Diff.Dy2];
        b       = bDiag - Cn*(Diff.Dy2*phi).^2;
    end

    A = [Ap Auv];
end                     