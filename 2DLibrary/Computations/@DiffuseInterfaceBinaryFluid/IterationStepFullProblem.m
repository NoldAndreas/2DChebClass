function IterationStepFullProblem(this,noIterations)
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
        
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;    
	nParticles     = this.optsPhys.nParticles;
    
    Diff           = this.IC.Diff;
    M              = this.IC.M;  
    Ind            = this.IC.Ind;
    IntSubArea     = this.IntSubArea;    
    y2Max          = this.optsNum.PhysArea.y2Max;
    

    Z    = zeros(M);
    IBB  = repmat(Ind.bound,2,1);  
    EYM  = eye(M);  EYMM  = eye(2*M);

    F       = false(M,1);   T       = true(M,1);    
    
    opts = struct('optsNum',this.optsNum,...
                  'optsPhys',this.optsPhys,...
                  'Comments',this.configName);
    
    [res,~,Parameters] = DataStorage([],@SolveBinaryFluid,...
                    opts,struct('initialGuess',GetInitialCondition(this)));
	   
    this.uv       = res.uv;
    this.p        = res.p;
    this.mu       = res.G;
    this.phi      = res.phi;    
    this.filename = Parameters.Filename;
    this.errors.errorIterations = res.errorHistory;
    
    %     this.errors.errorIterations = res.errorIterations;
    
    
    function res = SolveBinaryFluid(conf,in)
        [vec,errHistory] = NewtonMethod(in.initialGuess,@f,1e-6,noIterations,0.8);    
    
        res.uv  = vec([T;T;F;F;F]);
        res.phi = vec([F;F;T;F;F]);
        res.G   = vec([F;F;F;T;F]);
        res.p   = vec([F;F;F;F;T]);
        res.errorHistory = errHistory;
    end    
    function [v,A] = f(z)
        
        %[uv;phi;G;p] 
        uv  = z([T;T;F;F;F]);
        phi = z([F;F;T;F;F]);
        G   = z([F;F;F;T;F]);
        p   = z([F;F;F;F;T]);        
        
        phiM        = phi(Ind.left & Ind.top);
        phiP        = phi(Ind.right & Ind.top);  
        
        pM          = p(Ind.left & Ind.top);
        pP          = p(Ind.right & Ind.top);
        
        IntPathUpLow     = this.IC.borderTop.IntSc - this.IC.borderBottom.IntSc; 
        [fWP,fW,fWPP]    = DoublewellPotential(phi,Cn);
    

        % Continuity   %[uv;phi;G;p]
        A_cont         = [Diff.div,Z,Z,Z];    
        v_cont         = Diff.div*uv;

        % Momentum     %[uv;phi;G;p]         
         A_mom          = [Diff.LapVec,...
                           diag(repmat(fWP-Cn*Diff.Lap*phi,2,1))*Diff.grad/Cak + ...
                           diag(Diff.grad*phi)*repmat(diag(fWPP)-Cn*Diff.Lap,2,1)/Cak,...          
                           [Z;Z],...
                           - Diff.grad];
         v_mom          = Diff.LapVec*uv ...
                          + repmat(fWP-Cn*Diff.Lap*phi,2,1).*(Diff.grad*phi)/Cak ...
                          - Diff.grad*p;
        
        
        [v_G,A_G]   = Phasefield(this,uv,phi,G);
        [v_mu,A_mu] = ChemicalPotential(this,uv,phi,G);
                
       %% Boundary conditions [uv;phi;G;p]
        
       % (BC1) p = 0  at y1 = +/- infinity
       A_cont(Ind.left,:)           = 0;
       A_cont(Ind.left,[F;F;F;F;T]) = EYM(Ind.left,:);
       v_cont(Ind.left)             = p(Ind.left);
        
        % (BC2) [uv;phi;G;p]
       
       [fpW_mInf,fW_mInf] = DoublewellPotential(phiM,Cn);
       [fpW_pInf,fW_pInf] = DoublewellPotential(phiP,Cn);
     
       indBC2 = Ind.right & Ind.top;
     
       A_cont(Ind.right,:)           = 0;
       A_cont(Ind.right,[F;F;F;F;T]) = Diff.Dy2(Ind.right,:);
       v_cont(Ind.right,:)           = Diff.Dy2(Ind.right,:)*p;
       
       A_cont(indBC2,[T;T;F;F;F]) = IntPathUpLow*[Diff.Dy2 , Diff.Dy1];              
       
       A_cont(indBC2,[F;F;T;F;F])                 = -Cn/Cak*IntPathUpLow*(diag(Diff.Dy1*phi)*Diff.Dy2 + diag(Diff.Dy2*phi)*Diff.Dy1);
       A_cont(indBC2,[F;F;Ind.top&Ind.right;F;F]) = ...
                        A_cont(indBC2,[F;F;Ind.top&Ind.right;F;F]) + ...
                        y2Max*fpW_pInf/Cak;
       A_cont(indBC2,[F;F;Ind.top&Ind.left;F;F])  = ...
                        A_cont(indBC2,[F;F;Ind.top&Ind.left;F;F]) - ...
                        y2Max*fpW_mInf/Cak;
                    
       A_cont(indBC2,[F;F;F;F;Ind.top&Ind.right]) = -y2Max;  
       A_cont(indBC2,[F;F;F;F;Ind.top&Ind.left])  = y2Max;  
       
                    
       v_cont(indBC2) = IntPathUpLow*([Diff.Dy2 , Diff.Dy1]*uv ...
                              - Cn/Cak*((Diff.Dy1*phi).*(Diff.Dy2*phi)))...
                        +y2Max*((-pP + fW_pInf/Cak) - (-pM + fW_mInf/Cak));
        
        % (BC3) uv = uv_BC    
        [uvBound,a]            = GetBoundaryCondition(this);%,theta,phi);   
        A_mom(IBB,:)           = 0;
        A_mom(IBB,[T;T;F;F;F]) = EYMM(IBB,:);
        v_mom(IBB)             = uv(IBB) - uvBound(IBB);                
        
                
         Jint                                   = IntSubArea;        
         A_G(Ind.left & Ind.bottom,:)           = 0;
         A_G(Ind.left & Ind.bottom,[F;F;T;F;F]) = Jint;
         v_G(Ind.left & Ind.bottom)             = IntSubArea*phi - nParticles;                                 
        
        
        A = [A_mom;A_mu;A_G;A_cont];
        v = [v_mom;v_mu;v_G;v_cont];    
        
        DisplayError(v);
    end    
    function DisplayError(error)        
        PrintErrorPos(error([F;F;F;F;T]),'continuity equation',this.IC.Pts);
        PrintErrorPos(error([T;F;F;F;F]),'y1-momentum equation',this.IC.Pts);
        PrintErrorPos(error([F;T;F;F;F]),'y2-momentum equation',this.IC.Pts);                            
        PrintErrorPos(error([F;F;F;T;F]),'G equation',this.IC.Pts);
        PrintErrorPos(error([F;F;T;F;F]),'mu equation',this.IC.Pts);
    end

end