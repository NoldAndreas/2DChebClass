function [A,b] = FullStressTensorIJ(this,rho,i,j)
    % get matrices for
    % T = Cak*( eta*(grad(u) + grad(u)^T) + (zeta - 2/3*eta) div(u)*I ) +...
    %           + (W(rho) + Cn/2*|grad(rho)|^2 - mu*(rho+rho_m))*I - Cn*(grad(rho) X grad(rho))
    %   = A*[mu;uv] + b

    Cn    = this.optsPhys.Cn;
    Cak   = this.optsPhys.Cak;
    eta   = this.optsPhys.eta;            
    zeta  = this.optsPhys.zeta;
    rho_m = this.optsPhys.rho_m;
    Diff  = this.IC.Diff;
    M     = this.IC.M;

    [~,W]   = DoublewellPotential(rho,Cn);
    bDiag   = W + Cn/2*((Diff.Dy1*rho).^2 + (Diff.Dy2*rho).^2); %CahnHilliardFreeEnergy(rho,Cn,Diff);    
    if(i == 1 && j == 1)
        Amu     = -diag(rho+rho_m);
        Auv     = Cak*(eta*[2*Diff.Dy1 , zeros(M)] + (zeta - 2/3*eta)*Diff.div);
        b       = bDiag - Cn*(Diff.Dy1*rho).*(Diff.Dy1*rho);
    elseif((i==2 && j == 1)  || (i==1 && j == 2))
        Amu     = zeros(M);
        Auv     = Cak*[Diff.Dy2 , Diff.Dy1]*eta;
        b       = - Cn*(Diff.Dy1*rho).*(Diff.Dy2*rho);
    elseif(i==2 && j == 2)
        Amu     = -diag(rho+rho_m);
        Auv     = Cak*(2*eta*[zeros(M) , Diff.Dy2] + (zeta - 2/3)*Diff.div);
        b       = bDiag - Cn*(Diff.Dy2*rho).*(Diff.Dy2*rho);
    end

    A = [Amu Auv];

end                     