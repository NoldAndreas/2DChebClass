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
    
    if(length(this.optsPhys.thetaEq) == 1)
        Seppecher   = true;
    else
        Seppecher   = false;
    end
    solveSquare = true;
        
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
    
    if(Seppecher)
        if(isempty(this.theta))
            in.initialGuess = [0;0;pi/2;in.initialGuess];
        else
            in.initialGuess = [this.a;...
                               this.deltaX;...
                               this.theta;...
                               in.initialGuess];
        end
    end
    
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
        
        solveSquare = true;
        [vec,errHistory] = NewtonMethod(in.initialGuess,@f,1e-6,noIterations,0.8);    
        
        solveSquare = false;
        optsFS = optimoptions(@fsolve,'Jacobian','on','Display','iter','Algorithm','levenberg-marquardt');
        vec = fsolve(@f,vec,optsFS);    
    
        
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
            a      = z(1);
            deltaX = z(2);
            theta  = z(3);
            disp(['[a,deltaX,theta] = ',num2str(a),' , ',num2str(deltaX),' , ',num2str(theta*180/pi),'.']);
        
            z      = z(4:end);
        end

        uv  = z([T;T;F;F;F]);
        phi = z([F;F;T;F;F]);
        G   = z([F;F;F;T;F]);
        p   = z([F;F;F;F;T]);               
            
        [v_cont,A_cont] = Continuity(this,uv,phi,G,p);
        [v_mom,A_mom]   = Momentum(this,uv,phi,G,p);
        [v_G,A_G]       = Phasefield(this,uv,phi,G);
        [v_mu,A_mu]     = ChemicalPotential(this,uv,phi,G);
       
        A_particles            = zeros(1,5*M);
        A_particles([F;F;T;F]) = IntSubArea;
        v_particles            = IntSubArea*phi - nParticles; 
        
        if(solveSquare)
            lbC         = Ind.right & Ind.bottom;                
            A_G(lbC,:)  = A_particles;
            v_G(lbC,:)  = v_particles;                           
        end

        if(Seppecher)    
            A_cont = [zeros(M,3),A_cont];
            A_mom  = [zeros(2*M,3),A_mom];
            A_G    = [zeros(M,3),A_G];
            A_mu   = [zeros(M,3),A_mu];
    
            [v_mom(IBB),A_mom(IBB,:),v_mu(Ind.top),A_mu(Ind.top,:)]...        
                    = GetSeppecherBoundaryConditions(this,uv,phi,a,deltaX,theta);               
    
            [v_SeppAdd,A_SeppAdd] = GetSeppecherConditions(this,uv,phi,G,a,deltaX,theta);

            A_particles = [0,0,0,A_particles];
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
    end    
    function DisplayError(error)
        
        if(Seppecher)
            PrintErrorPos(error(1),'consistent mass influx');
            PrintErrorPos(error(2),'zero density at interface');        
            PrintErrorPos(error(3),'chemical potential at -inf');
        
            error = error(4:end);
        end
        PrintErrorPos(error([F;F;F;F;T]),'continuity equation',this.IC.Pts);
        PrintErrorPos(error([T;F;F;F;F]),'y1-momentum equation',this.IC.Pts);
        PrintErrorPos(error([F;T;F;F;F]),'y2-momentum equation',this.IC.Pts);                            
        PrintErrorPos(error([F;F;F;T;F]),'G equation',this.IC.Pts);
        PrintErrorPos(error([F;F;T;F;F]),'mu equation',this.IC.Pts);
    end

end