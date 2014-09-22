function [A,b] = ContinuityMomentumEqs(this,phi)
    %Continuitiy: div(phi*uv) = 0
    %Momentum: - (phi+phi_m)*grad(mu) +
    %          ... +grad(phi)*(W'(phi) - mu - Cn*Lap(phi) )
    %          + Cak*(eta*Lap(uv) + (zeta + eta/3)*grad(div(uv)) )

    %
    % A*[mu;uv] = b corresponds to momentum and continuity Eq. for given phi               
    phi_m          = this.optsPhys.phi_m;
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;
    eta            = this.optsPhys.eta;
    zeta           = this.optsPhys.zeta;
    Diff           = this.IC.Diff;
    M              = this.IC.M;

    phi_f          = phi + phi_m;
    phi_f2         = repmat(phi_f,2,1);
    gradRho_T      = Diff.grad*phi_f;

    aT             = zeros(M,2*M);
    aT(:,1:M)      = diag(gradRho_T(1:M));
    aT(:,1+M:end)  = diag(gradRho_T(1+M:end));

    A_cont_mu      = zeros(M);
    %A_cont_uv      = Diff.div*diag(phi_f2);%aT + diag(phi_f)*Diff.div; 
    A_cont_uv      = aT + diag(phi_f)*Diff.div; 

    %A_mom_mu       = -diag(phi_f2)*Diff.grad - [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)];
    A_mom_mu       = - Diff.grad*diag(phi_f);
    A_mom_uv       = Cak*(eta*Diff.LapVec + (zeta + eta/3)*Diff.gradDiv);

    A_cont         = [A_cont_mu,A_cont_uv];
    A_mom          = [A_mom_mu, A_mom_uv];
    A              = [A_cont;A_mom];  

    b              = zeros(3*M,1);            
    ys             = DoublewellPotential(phi,Cn) - Cn*(Diff.Lap*phi);

    b(1+M:end)     = - repmat(ys,2,1).*(Diff.grad*phi); 
end        