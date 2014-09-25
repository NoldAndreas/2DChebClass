function IterationStepFullProblem(this,vec)
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
    s              = this.s;
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;    
    m              = this.optsPhys.mobility;

    Diff           = this.IC.Diff;
    M              = this.IC.M;  
    Ind            = this.IC.Ind;

    Z    = zeros(M); ZM = zeros(M,1);
    IBB  = repmat(Ind.bound,2,1);  
    EYM  = eye(M);  EYMM  = eye(2*M);

    F       = false(M,1);   T       = true(M,1);
	
    vec = NewtonMethod(vec,@f,1e-6,100,0.1);
    
    function [v,A] = f(z)
        
        %[uv;phi;G;p] 
        uv  = z([T;T;F;F;F]);
        phi = z([F;F;T;F;F]);
        G   = z([F;F;F;T;F]);
        p   = z([F;F;F;F;T]);        
   
        gradPhiM       = [diag(Diff.Dy1*phi);diag(Diff.Dy2*phi)];
        gradPhiMT      = [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)];

        % continuity
        A_cont         = [Diff.div,Z,Z,Z];    
        v_cont         = A_cont*vec;

        % momentum
        A_mom          = [Diff.LapVec,...
                          diag(repmat(G+s*phi,2,1))*Diff.grad/Cak+s*gradPhiM/Cak,...
                            gradPhiM/Cak,...
                            - Diff.grad];
        v_mom          = Diff.LapVec*uv ...
                         - repmat(G+s*phi,2,1).*(Diff.grad*phi)/Cak ...
                         - Diff.grad*p;

        % Phasefield
        A_G            = [-gradPhiMT,s*m*Diff.Lap-gradPhiMT*Diff.grad,m*Diff.Lap,Z];
        v_G            = m*Diff.Lap*(G + s*phi) - gradPhiMT*uv;    

        % Chemical Potential
        [fWP,~,fWPP]   = DoublewellPotential(phi,Cn);    
        A_mu           = [Z,Z,-s-Cn*Diff.Lap+diag(fWPP),-EYM,Z];
        v_mu           = fWP - s*phi - Cn*Diff.Lap*phi - G;

        %% Boundary conditions [uv;phi;G;p]

        % (BC1) uv = uv_BC    
        [uvBound,a]            = GetBoundaryCondition(this);%,theta,phi);   
        A_mom(IBB,:)           = 0;
        A_mom(IBB,[T;T;F;F;F]) = EYMM(IBB,:);
        v_mom(IBB)             = uv(IBB) - uvBound(IBB);

        % (BC2) nu*grad(phi) = 0
        A_mu(Ind.top|Ind.bottom,:)           = 0;
        A_mu(Ind.top|Ind.bottom,[F;F;T;F;F]) = Diff.Dy2(Ind.top|Ind.bottom,:);    
        v_mu(Ind.top|Ind.bottom)             = Diff.Dy2(Ind.top|Ind.bottom,:)*phi;

        % (BC3) nu*grad(G) = 0
        A_G(Ind.top|Ind.bottom,:)           = 0;
        A_G(Ind.top|Ind.bottom,[F;F;F;T;F]) = Diff.Dy2(Ind.top|Ind.bottom,:);    
        v_G(Ind.top|Ind.bottom)             = Diff.Dy2(Ind.top|Ind.bottom,:)*G;    

        % (BC4) p = 0  at y1 = +/- infinity
        A_cont(Ind.left|Ind.right,:)           = 0;
        A_cont(Ind.left|Ind.right,[F;F;F;F;T]) = EYM(Ind.left|Ind.right,:);
        v_cont(Ind.left|Ind.right)             = p(Ind.left|Ind.right);

        A = [A_mom;A_mu;A_G;A_cont];
        v = [v_mom;v_mu;v_G;v_cont];
    
    end
    
end