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
    zeta           = this.optsPhys.zeta;
    phi_m          = this.optsPhys.phi_m;            
    
    Diff           = this.IC.Diff;
    M              = this.IC.M;  
    Ind            = this.IC.Ind;    
    tauM           = Cak*(Diff.LapVec + (zeta + 1/3)*Diff.gradDiv);
    
    lbC = Ind.left & Ind.bottom;
    rbC = Ind.right & Ind.bottom;

    Z    = zeros(M);
    IBB  = repmat(Ind.bound,2,1);  
    EYMM  = eye(2*M);

    F       = false(M,1);   T  = true(M,1);    
    
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
    
        [v_cont,A_cont]  = Continuity(this,uv,phi,G);               
        [v_mom,A_mom]    = Momentum(this,uv,phi,G);      
        [v_mu,A_mu]      = ChemicalPotential(this,uv,phi,G);
                        
        
       %% Boundary conditions [uv;phi;G]                              
        
        % (BC3) uv = uv_BC    
        [uvBound,a]            = GetBoundaryCondition(this);%,theta,phi);   
        A_mom(IBB,:)           = 0;
        A_mom(IBB,[T;T;F;F])   = EYMM(IBB,:);
        v_mom(IBB)             = uv(IBB) - uvBound(IBB);
        
        % (BC4) nu*grad(phi) = 0
        A_cont(Ind.top,:)           = 0;
        A_cont(Ind.top,[F;F;T;F;F]) = Diff.Dy2(Ind.top,:);    
        v_cont(Ind.top)             = Diff.Dy2(Ind.top,:)*phi;
        
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