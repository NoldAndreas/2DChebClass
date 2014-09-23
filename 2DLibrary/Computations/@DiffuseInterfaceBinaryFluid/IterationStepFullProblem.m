function IterationStepFullProblem(this)
    %Continuitiy: 0 = div(uv)
    %Momentum:    0 = -Grad(p) + (G+s*phi)*grad(phi)/Cak + Lap(uv)
    %Phasefield   0 = m*Lap(G + s*phi) - u*grad(phi)
    %ChemPot      0 = f_w' - s*phi - Cn*Lap(phi) - G
    %
    % 
    % (BC1) uv = uv_BC    
    % (BC2) nu*grad(phi) = 0
    % (BC3) nu*grad(G) = 0
    % (BC4) p = 0  at y1 = +/- infinity
    
    % A*[uv;phi;G;p] = b corresponds to momentum and continuity Eq. for given phasefield phi               

    s = 0.5;
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;    
    m              = this.optsPhys.mobility;
    
    Diff           = this.IC.Diff;
    M              = this.IC.M;            

    Z    = zeros(M); ZM = zeros(M,1);
    IBB  = repmat(Ind.bound,2,1);  
    EYM  = eye(M);  EYMM  = eye(2*M);
    
    
    mu = this.mu; uv = this.uv; phi = this.phi; p = this.p;
    G  = mu - s*phi;
    vec = [uv;phi;mu;p];
    
    gradPhiM       = [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)];
    gradPhiMT      = [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)];
      
    % continuity
    A_cont         = [Diff.div,Z,Z,Z];    
    f_cont         = A_cont*vec;
    
    % momentum
    A_mom          = [Diff.LapVec,...
                      diag(G+s*phi)*Diff.grad/Cak+s*gradPhiM/Cak,...
                        gradPhiM/Cak,...
                        - Diff.grad];
    f_mom          = Diff.LapVec*uv ...
                     - repmat(G+s*phi,2,1).*(Diff.grad*phi)/Cak ...
                     - Diff.grad*p;
    
    % Phasefield
    A_G            = [-gradPhiMT,s*m*Diff.Lap-gradPhiMT*Diff.grad,m*Diff.Lap,Z];
    f_G            = m*Diff.Lap*(G + s*phi) - gradPhiMT*uv;    
    
    % Chemical Potential
    [fWP,~,fWPP]   = DoubleWellPotential(phi,Cn);    
    A_mu           = [Z,Z,-s-Cn*Diff.Lap+diag(fWPP),-EYM,Z];
    f_mu           = fWP - s*phi - Cn*Diff.Lap*phi - G;
    
    %% Boundary conditions
    
    % BC1
    [uvBound,a]  = GetBoundaryCondition(this,theta,phi);   
    A_mom(IBB,:) = EYMM(IBB,:)
    f_mom(IBB)   = uv(IBB) - uvBound(IBB);
	
    % BC2    
    A_mu(Ind.top|Ind.bottom,:) = Diff.Dy2(Ind.top|Ind.bottom,:);    
    f_mu(Ind.top|Ind.bottom)   = Diff.Dy2(Ind.top|Ind.bottom,:)*phi;
    
    
end