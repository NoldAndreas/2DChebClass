function [A,b] = ContinuityMomentumEqs_mu_p_uv(this,phi)
    %Continuitiy: 0 = div(uv)
    %Momentum:    0 = -Grad(p) + Lap(uv) + 1/Cak * mu*grad(phi)
    %Chem Pot:    0 = u*grad(phi) - m*Lap(mu)
    %
    % A*[p;mu;uv] = b corresponds to momentum and continuity Eq. for given phasefield phi               
                  
    Cak            = this.optsPhys.Cak;    
    m              = this.optsPhys.mobility;
    
    Diff           = this.IC.Diff;
    M              = this.IC.M;            
    
    gradPhiT       = [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)];
      
    A_cont         = [zeros(M),zeros(M),Diff.div];        
    A_mom          = [- Diff.grad,gradPhiT/Cak,Diff.LapVec];
    A_mu           = [zeros(M),-m*Diff.Lap,gradPhiT];
    
    
    A              = [A_cont;A_mom;A_mu];  
    b              = zeros(4*M,1);                
end        