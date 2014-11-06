function [A,b] = FullStressTensorIJ(this,phi,i,j)
    % get matrices for
    % T = Cak*( eta*(grad(u) + grad(u)^T) + (zeta - 2/3*eta) div(u)*I ) +...
    %           + (W(phi) + Cn/2*|grad(phi)|^2 - mu*(phi+phi_m))*I - Cn*(grad(phi) X grad(phi))
    %   = A*[mu;uv] + b

    Cn    = this.optsPhys.Cn;
    Cak   = this.optsPhys.Cak;               
    zeta  = this.optsPhys.zeta;
    phi_m = this.optsPhys.phi_m;
    Diff  = this.IDC.Diff;
    M     = this.IDC.M;

    [~,W]   = DoublewellPotential(phi,Cn);
    bDiag   = W + Cn/2*((Diff.Dy1*phi).^2 + (Diff.Dy2*phi).^2); %CahnHilliardFreeEnergy(phi,Cn,Diff);    
    if(i == 1 && j == 1)
        Amu     = -diag(phi+phi_m);
        Auv     = Cak*([2*Diff.Dy1 , zeros(M)] + (zeta - 2/3)*Diff.div);
        b       = bDiag - Cn*(Diff.Dy1*phi).*(Diff.Dy1*phi);
    elseif((i==2 && j == 1)  || (i==1 && j == 2))
        Amu     = zeros(M);
        Auv     = Cak*[Diff.Dy2 , Diff.Dy1];
        b       = - Cn*(Diff.Dy1*phi).*(Diff.Dy2*phi);
    elseif(i==2 && j == 2)
        Amu     = -diag(phi+phi_m);
        Auv     = Cak*(2*[zeros(M) , Diff.Dy2] + (zeta - 2/3)*Diff.div);
        b       = bDiag - Cn*(Diff.Dy2*phi).*(Diff.Dy2*phi);
    end

    A = [Amu Auv];

end                     