function IterationStepFullProblem_old(this,noIterations)
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
        
    solveSquare = true;
    
    nParticles     = this.optsPhys.nParticles;
        
    M              = this.IC.M;  
    Ind            = this.IC.Ind;
    IntSubArea     = this.IntSubArea;        
        
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
        
        solveSquare = true;
        [vec,errHistory] = NewtonMethod(in.initialGuess,@f,1e-6,noIterations,0.8);    
        %vec = in.initialGuess;  errHistory = [];
        
        solveSquare = false;
        optsFS = optimoptions(@fsolve,'Jacobian','on','Display','iter','Algorithm','levenberg-marquardt');
        vec = fsolve(@f,vec,optsFS);
        
        %[vec,errHistory] = NewtonMethod(in.initialGuess,@f,1e-6,noIterations,0.8);    
    
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

        [v_cont,A_cont] = Continuity(this,uv,phi,G,p);
        [v_mom,A_mom]   = Momentum(this,uv,phi,G,p);
        [v_G,A_G]       = Phasefield(this,uv,phi,G);
        [v_mu,A_mu]     = ChemicalPotential(this,phi,G);
                
       %% Boundary conditions [uv;phi;G;p]
        

        % (BC3) uv = uv_BC    
        [v_mom(IBB),A_mom(IBB,:)] = GetVelBC(this,uv);
%         [uvBound,a]            = GetBoundaryCondition(this);%,theta,phi);   
%         A_mom(IBB,:)           = 0;
%         A_mom(IBB,[T;T;F;F;F]) = EYMM(IBB,:);
%         v_mom(IBB)             = uv(IBB) - uvBound(IBB);                
        
        A_particles            = zeros(1,5*M);
        A_particles([F;F;T;F;F]) = IntSubArea;
        v_particles              = IntSubArea*phi - nParticles; 
        if(solveSquare)                    
            A_G(Ind.right & Ind.bottom,:)   = A_particles;            
            v_G(Ind.right & Ind.bottom)     = v_particles;   
            
            A = [A_mom;A_mu;A_G;A_cont];
            v = [v_mom;v_mu;v_G;v_cont];            
        else
            A = [A_mom;A_mu;A_G;A_cont;A_particles];
            v = [v_mom;v_mu;v_G;v_cont;v_particles];    
        end
                        
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