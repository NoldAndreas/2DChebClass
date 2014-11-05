function IterationStepFullProblem(this,opts)
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
    if((nargin == 1) || ~isfield(opts,'noIterations'))
        noIterations = 20;
    else
        noIterations = opts.noIterations;
    end
    
    if((nargin == 1) || ~isfield(opts,'lambda'))
        lambda = 0.8;
    else
        lambda = opts.lambda;
    end
    
	if((nargin == 1) || ~isfield(opts,'solveSquare'))
        solveSquare = true;
    else
        solveSquare = opts.solveSquare;
    end
    
    Seppecher   = IsSeppecher(this);    
    %solveSquare = true;
        
    nParticles     = this.optsPhys.nParticles;            
    M              = this.IC.M;  
    Ind            = this.IC.Ind;
    IntSubArea     = this.IntSubArea;        
    
    
    IBB  = repmat(Ind.bound,2,1);      
    F       = false(M,1);   T       = true(M,1);        
    
    opts = struct('optsNum',this.optsNum,...
                  'optsPhys',this.optsPhys,...
                  'Comments',this.configName);
              
    
	in = struct('initialGuess',GetInitialCondition(this));
        
    [res,~,Parameters] = DataStorage([],@SolveSingleFluid,...
                    opts,in);
	   
    this.uv       = res.uv;
    this.mu       = res.mu;
    this.phi      = res.phi;
    this.p        = res.p;
    if(Seppecher)
        this.a        = res.a;        
        this.deltaX   = res.deltaX;
        this.theta    = res.theta;
    end
    
    this.filename = Parameters.Filename;
    this.errors.errorIterations = res.errorHistory;        
    
    
    function res = SolveSingleFluid(conf,in)
        
        vec = in.initialGuess;
        %solveSquare = true;
        if(solveSquare)
            [vec,errHistory] = NewtonMethod(vec,@f,1e-6,noIterations,lambda);    
        else        
            %solveSquare = false;
            optsFS = optimoptions(@fsolve,'Jacobian','on','Display','iter','Algorithm','levenberg-marquardt');
            vec = fsolve(@f,vec,optsFS);    
            errHistory = [];
        end
    
        
        %[uv;phi;mu] 
        if(Seppecher)
            res.a      = vec(1);
            res.deltaX = vec(2);
            res.theta  = vec(3);
            vec        = vec(4:end);
        end
        res.uv     = vec([T;T;F;F;F]);
        res.phi    = vec([F;F;T;F;F]);
        res.mu     = vec([F;F;F;T;F]);
        res.p      = vec([F;F;F;F;T]);
        
        res.errorHistory = errHistory;
    end    
    function [v,A] = f(z)
        
        %[uv;phi;G;p]       
        if(Seppecher)
            a      = [z(1),0];
            deltaX = z(2);
            theta  = z(3);
            disp(['[a,deltaX,theta] = ',num2str(a(1)),' , ',num2str(a(2)),' , ',num2str(deltaX),' , ',num2str(theta*180/pi),'.']);
        
            z      = z(4:end);
        else
            a = []; deltaX = []; theta = [];                    
        end

        uv  = z([T;T;F;F;F]);
        phi = z([F;F;T;F;F]);
        G   = z([F;F;F;T;F]);
        p   = z([F;F;F;F;T]);               
            
        [v_cont,A_cont] = Continuity(this,uv,phi,G,p);
        [v_mom,A_mom]   = Momentum(this,uv,phi,G,p);
        [v_G,A_G]       = Phasefield(this,uv,phi,G);
        [v_mu,A_mu]     = ChemicalPotential(this,phi,G);
       
        A_particles              = zeros(1,5*M);
        A_particles([F;F;T;F;F]) = IntSubArea;
        v_particles              = IntSubArea*phi - nParticles;         

        if(Seppecher)    
            A_cont      = [zeros(M,4),A_cont];
            A_mom       = [zeros(2*M,4),A_mom];
            A_G         = [zeros(M,4),A_G];
            A_mu        = [zeros(M,4),A_mu];
            A_particles = [0,0,0,0,A_particles];
            
            [v_SeppAdd,A_SeppAdd] = GetSeppecherConditions(this,uv,phi,G,p,deltaX,theta);
        end
        
        [v_mu(Ind.top|Ind.bottom),A_mu(Ind.top|Ind.bottom,:)] = GetPhiBC(this,phi,theta);          
        [v_mom(IBB),A_mom(IBB,:)]                             = GetVelBC(this,uv,a,deltaX,theta);                    
        [v_G(Ind.bound),A_G(Ind.bound,:)]                     = GetChemPotBC(this,uv,phi,G,theta);
        
        if(solveSquare)
            rbC         = Ind.right & Ind.bottom;                
            A_G(rbC,:)  = A_particles;
            v_G(rbC,:)  = v_particles;                           
        end
        
        A = [A_mom;A_mu;A_G;A_cont];
        v = [v_mom;v_mu;v_G;v_cont];

        if(Seppecher)
            A  = [A_SeppAdd;A];
            v  = [v_SeppAdd;v];
        end       

        if(~solveSquare)
            A = [A;A_particles];
            v = [v;v_particles];                           
        end
        
        DisplayError(v);
                
        A = A(2:end,2:end);
        v = v(2:end);
        
        %A = A([1,3:end],[1,3:end]);
        %v = v([1,3:end]);
    end    
    function DisplayError(error)
        
        if(Seppecher)
            PrintErrorPos(error(1),'consistent phase mass');
            PrintErrorPos(error(2),'consistent mass');
            PrintErrorPos(error(3),'zero density at interface');        
            PrintErrorPos(error(4),'chemical potential at -inf');
        
            error = error(5:end);
        end
        PrintErrorPos(error([F;F;F;F;T]),'continuity equation',this.IC.Pts);
        PrintErrorPos(error([T;F;F;F;F]),'y1-momentum equation',this.IC.Pts);
        PrintErrorPos(error([F;T;F;F;F]),'y2-momentum equation',this.IC.Pts);                            
        PrintErrorPos(error([F;F;F;T;F]),'G equation',this.IC.Pts);
        PrintErrorPos(error([F;F;T;F;F]),'mu equation',this.IC.Pts);
    end

end