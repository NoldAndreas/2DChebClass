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
    if(nargin == 1)
        noIterations = 20;
    end
    
    Cn             = this.optsPhys.Cn;
    Cak            = this.optsPhys.Cak;
    eta            = this.optsPhys.eta;
    zeta           = this.optsPhys.zeta;
    phi_m          = this.optsPhys.phi_m;    
    nParticles     = this.optsPhys.nParticles;
    
    Diff           = this.IC.Diff;
    M              = this.IC.M;  
    Ind            = this.IC.Ind;
    IntSubArea     = this.IntSubArea;    
    y2Max          = this.optsNum.PhysArea.y2Max;
    
    IntPathUpLow   = this.IC.borderTop.IntSc - this.IC.borderBottom.IntSc; 
    tauM           = Cak*(eta*Diff.LapVec + (zeta + eta/3)*Diff.gradDiv);
    
    lbC = Ind.left & Ind.bottom;
    rbC = Ind.right & Ind.bottom;

    Z    = zeros(M);
    IBB  = repmat(Ind.bound,2,1);  
    EYM  = eye(M);  EYMM  = eye(2*M);

    F       = false(M,1);   T       = true(M,1);    
    
    opts = struct('optsNum',this.optsNum,...
                  'optsPhys',this.optsPhys,...
                  'Comments',this.configName);
    
    [res,~,Parameters] = DataStorage([],@SolveSingleFluid,...
                    opts,struct('initialGuess',GetInitialCondition(this)));
	   
    this.uv       = res.uv;
    this.mu       = res.mu;
    this.phi      = res.phi;
    this.filename = Parameters.Filename;
    this.errors.errorIterations = res.errorHistory;
    
    %     this.errors.errorIterations = res.errorIterations;
    
    
    function res = SolveSingleFluid(conf,in)
        [vec,errHistory] = NewtonMethod(in.initialGuess,@f,1e-6,noIterations,0.8);    
    
        %[uv;phi;mu] 
        res.uv  = vec([T;T;F;F;F]);
        res.phi = vec([F;F;T;F;F]);
        res.mu  = vec([F;F;F;T;F]);
        
        res.errorHistory = errHistory;
    end    
    function [v,A] = f(z)
        
        %[uv;phi;G;p] 
        uv  = z([T;T;F;F]);
        phi = z([F;F;T;F]);
        G   = z([F;F;F;T]);        
   
        uvTdiag        = [diag(uv(1:end/2)),diag(uv(end/2+1:end))];        
        [fWP,fW,fWPP]    = DoublewellPotential(phi,Cn);
    

        % Continuity   %[uv;phi;G]
        A_cont         = [diag(phi+phi_m)*Diff.div + [diag(Diff.Dy1*phi),diag(Diff.Dy2*phi)]...
                          diag(Diff.div*uv)+uvTdiag*Diff.grad,...
                          Z];
        v_cont         = Diff.div*(uv.*repmat(phi+phi_m,2,1));

       % Momentum     %[uv;phi;G]       
       A_mom          = [tauM,...
                         -[diag(Diff.Dy1*G);diag(Diff.Dy2*G)],...
                         - diag(repmat(phi+phi_m,2,1))*Diff.grad];

       v_mom         = tauM*uv - repmat(phi+phi_m,2,1).*(Diff.grad*G);
        
       % Chemical Potential %[uv;phi;G]  
       A_mu           = [Z,Z,...
                         diag(fWPP)-Cn*Diff.Lap,...
                         -EYM];
       v_mu           = fWP - Cn*Diff.Lap*phi - G;
        
       %% Boundary conditions [uv;phi;G]
        
       % (BC1) p = 0  at y1 = +/- infinity
       A_cont(Ind.left|Ind.right,:)         = 0;
       A_cont(Ind.left|Ind.right,[F;F;F;T]) = Diff.Dy2(Ind.left|Ind.right,:);
       v_cont(Ind.left|Ind.right,:)         = Diff.Dy2(Ind.left|Ind.right,:)*G;
                               
       A_cont(lbC,:)         = 0;
       A_cont(lbC,[F;F;T;F]) = IntSubArea;
       v_cont(lbC)           = IntSubArea*phi - nParticles;                 
               
       % (BC2) [uv;phi;G]                   
       A_cont(rbC,[T;T;F;F])   = IntPathUpLow*[Diff.Dy2 , Diff.Dy1]*eta;       
       A_cont(rbC,[F;F;T;F])   = -Cn/Cak*IntPathUpLow*(diag(Diff.Dy1*phi)*Diff.Dy2 + diag(Diff.Dy2*phi)*Diff.Dy1);
       
       A_cont(rbC,[F;F;rbC;F]) = A_cont(rbC,[F;F;rbC;F]) + y2Max*(G(rbC) + fWP(rbC))/Cak;
       A_cont(rbC,[F;F;lbC;F]) = A_cont(rbC,[F;F;lbC;F]) - y2Max*(G(lbC) + fWP(lbC))/Cak;       
       A_cont(rbC,[F;F;F;rbC]) = A_cont(rbC,[F;F;F;rbC]) + y2Max*(phi(rbC) + phi_m)/Cak;
       A_cont(rbC,[F;F;F;lbC]) = A_cont(rbC,[F;F;F;lbC]) - y2Max*(phi(lbC) + phi_m)/Cak;
                    
       v_cont(rbC) = IntPathUpLow*(eta*[Diff.Dy2 , Diff.Dy1]*uv ...
                              - Cn/Cak*((Diff.Dy1*phi).*(Diff.Dy2*phi)))...
                        +y2Max*(((phi(rbC)+phi_m)*G(rbC) + fW(rbC)) - ...
                                ((phi(lbC)+phi_m)*G(lbC) + fW(lbC)))/Cak;
        
        % (BC3) uv = uv_BC    
        [uvBound,a]            = GetBoundaryCondition(this);%,theta,phi);   
        A_mom(IBB,:)           = 0;
        A_mom(IBB,[T;T;F;F;F]) = EYMM(IBB,:);
        v_mom(IBB)             = uv(IBB) - uvBound(IBB);
        
        % (BC4) nu*grad(phi) = 0
        A_mu(Ind.top|Ind.bottom,:)           = 0;
        A_mu(Ind.top|Ind.bottom,[F;F;T;F;F]) = Diff.Dy2(Ind.top|Ind.bottom,:);    
        v_mu(Ind.top|Ind.bottom)             = Diff.Dy2(Ind.top|Ind.bottom,:)*phi;
                                                
        
        A = [A_mom;A_mu;A_cont];
        v = [v_mom;v_mu;v_cont];    
        
        DisplayError(v);
    end    
    function DisplayError(error)        
        PrintErrorPos(error([F;F;F;T]),'continuity equation',this.IC.Pts);
        PrintErrorPos(error([T;F;F;F]),'y1-momentum equation',this.IC.Pts);
        PrintErrorPos(error([F;T;F;F]),'y2-momentum equation',this.IC.Pts);                                    
        PrintErrorPos(error([F;F;T;F]),'mu equation',this.IC.Pts);
    end

end