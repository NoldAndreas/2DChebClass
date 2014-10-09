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
    
    s              = this.s;
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;    
    m              = this.optsPhys.mobility;
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
    this.mu       = res.G + s*res.phi;
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
   
        uvTdiag        = [diag(uv(1:end/2)),diag(uv(end/2+1:end))];        
        gradPhiMT      = [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)];
        
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
        
        % Phasefield   %[uv;phi;G;p]        
        A_G            = [-gradPhiMT,s*m*Diff.Lap-uvTdiag*Diff.grad,m*Diff.Lap,Z];
        v_G            = m*Diff.Lap*(G + s*phi) - gradPhiMT*uv;    

        % Chemical Potential        
        A_mu           = [Z,Z,diag(fWPP)-Cn*Diff.Lap-s*EYM,-EYM,Z];
        v_mu           = fWP - Cn*Diff.Lap*phi - (G + s*phi);
        
        lr = Ind.left|Ind.right;
        A_mu(lr,:) = 0;
        fwD        = diag(fWPP);
        
        A_mu(lr,:) = [Z(lr,:),...
                      Z(lr,:),...
                      fwD(lr,:)-Cn*Diff.DDy2(lr,:)-s*EYM(lr,:),...
                      -EYM(lr,:),...
                      Z(lr,:)];
        v_mu(lr)   = fWP(lr) - Cn*Diff.DDy2(lr,:)*phi ...
                             - (G(lr) + s*phi(lr));

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
        
        % (BC4) nu*grad(G) = 0
        A_G(Ind.top|Ind.bottom,:)           = 0;
        A_G(Ind.top|Ind.bottom,[F;F;F;T;F]) = Diff.Dy2(Ind.top|Ind.bottom,:);    
        v_G(Ind.top|Ind.bottom)             = Diff.Dy2(Ind.top|Ind.bottom,:)*G;                                    
        
        % (BC5) 0 = f_w' - s*phi - Cn*Lap(phi) - G at y1 = +/- infintity
        ind_LR         = Ind.left | Ind.right;
        fwPP_D         = diag(fWPP);
        A_G(ind_LR,:)  = [Z(ind_LR,:),Z(ind_LR,:),...
                          fwPP_D(ind_LR,:)-s*EYM(ind_LR,:),...%-Cn*Diff.Lap(ind_LR,:),...
                          -EYM(ind_LR,:),...
                          Z(ind_LR,:)];                
        v_G(ind_LR)    = fWP(ind_LR) - Cn*Diff.Lap(ind_LR,:)*phi ...
                                    - (G(ind_LR) + s*phi(ind_LR));

        % (BC6) nu*grad(phi) = 0
        A_mu(Ind.top|Ind.bottom,:)           = 0;
        A_mu(Ind.top|Ind.bottom,[F;F;T;F;F]) = Diff.Dy2(Ind.top|Ind.bottom,:);    
        v_mu(Ind.top|Ind.bottom)             = Diff.Dy2(Ind.top|Ind.bottom,:)*phi;
                                
        
        %(BC7) Mass in subarea:        
         Jint                      = IntSubArea;        
         A_mu(Ind.left & Ind.top,:)           = 0;
         A_mu(Ind.left & Ind.top,[F;F;T;F;F]) = Jint;
         v_mu(Ind.left & Ind.top)  = IntSubArea*phi - nParticles;                 
        
         %(BC8) [uv;phi;G;p]
        A_mu(Ind.right,:)           = 0;
        A_mu(Ind.right,[F;F;T;F;F]) = Diff.Dy2(Ind.right,:);
        v_mu(Ind.right,:)           = Diff.Dy2(Ind.right,:)*phi;
        
        hM  = - m*this.IC.borderBottom.IntSc*Diff.DDy2;        
        A_mu(Ind.right & Ind.bottom,:)           = 0;
        A_mu(Ind.right & Ind.bottom,[F;T;F;F;F]) = this.IC.borderBottom.IntSc*(diag(phi)*Diff.Dy2);
        A_mu(Ind.right & Ind.bottom,[F;F;T;F;F]) = + EYM(Ind.right & Ind.bottom,:) ...
                                                - EYM(Ind.left & Ind.bottom,:) ...
                                              + s*hM...
                                              + this.IC.borderBottom.IntSc*(Diff.Dy2*(uv(1+end/2:end)));
        A_mu(Ind.right & Ind.bottom,[F;F;F;T;F]) = hM;
        v_mu(Ind.right & Ind.bottom)             = phi(Ind.right & Ind.bottom) ...
                                               - phi(Ind.left & Ind.bottom) ...
                                               + this.IC.borderBottom.IntSc*(phi.*(Diff.Dy2*(uv(1+end/2:end))))...
                                               + hM*(G + s*phi);
        
        
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