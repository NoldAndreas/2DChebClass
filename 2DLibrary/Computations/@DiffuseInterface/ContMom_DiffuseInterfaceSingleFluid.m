function [A,b] = ContMom_DiffuseInterfaceSingleFluid(this,rho)
    %Continuitiy: div(rho*uv) = 0
    %Momentum: - (rho+rho_m)*grad(mu) +
    %          ... +grad(rho)*(W'(rho) - mu - Cn*Lap(rho) )
    %          + Cak*(eta*Lap(uv) + (zeta + eta/3)*grad(div(uv)) ) +...

    %
    % A*[mu;uv] = b corresponds to momentum and continuity Eq. for given rho               
    rho_m          = this.optsPhys.rho_m;
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;
    eta            = this.optsPhys.eta;
    zeta           = this.optsPhys.zeta;
    Diff           = this.IC.Diff;
    M              = this.IC.M;

    rho_f          = rho + rho_m;
    rho_f2         = repmat(rho_f,2,1);
    gradRho_T      = Diff.grad*rho_f;

    aT             = zeros(M,2*M);
    aT(:,1:M)      = diag(gradRho_T(1:M));
    aT(:,1+M:end)  = diag(gradRho_T(1+M:end));

    A_cont_mu      = zeros(M);
    A_cont_uv      = aT + diag(rho_f)*Diff.div; 

    %A_mom_mu       = -diag(rho_f2)*Diff.grad - [diag(Diff.Dy1*rho);diag(Diff.Dy2*rho)];
    A_mom_mu       = - Diff.grad*diag(rho_f);
    A_mom_uv       = Cak*(eta*Diff.LapVec + (zeta + eta/3)*Diff.gradDiv);

    A_cont         = [A_cont_mu,A_cont_uv];
    A_mom          = [A_mom_mu, A_mom_uv];
    A              = [A_cont;A_mom];  

    b              = zeros(3*M,1);            
    ys             = DoublewellPotential(rho,Cn) - Cn*(Diff.Lap*rho);

    b(1+M:end)     = - repmat(ys,2,1).*(Diff.grad*rho); 
end        