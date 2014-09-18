function [A,b] = ContMom_DiffuseInterfaceBinaryFluid(this,phi)
    %Continuitiy: div(rho*uv) = 0
    %Momentum: - (rho+rho_m)*grad(mu) +
    %          ... +grad(rho)*(W'(rho) - mu - Cn*Lap(rho) )
    %          + Cak*(eta*Lap(uv) + (zeta + eta/3)*grad(div(uv)) )

    %
    % A*[p;uv] = b corresponds to momentum and continuity Eq. for given phasefield phi               
    
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;    
    Diff           = this.IC.Diff;
    M              = this.IC.M;            

    A_cont_p       = zeros(M);    
    A_cont_uv      = Diff.div; 
    
    A_mom_p        = - Diff.grad;
    A_mom_uv       = Cak*Diff.LapVec;

    A_cont         = [A_cont_p,A_cont_uv];
    A_mom          = [A_mom_p, A_mom_uv];
    A              = [A_cont;A_mom];  

    b              = zeros(3*M,1);            
    ys             = DoublewellPotential(phi,Cn) - Cn*(Diff.Lap*rho);
    b(1+M:end)     = - repmat(ys,2,1).*(Diff.grad*phi); 
end        