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
    
    if(isfield(this.optsNum,'optsSolv'))
        opts = this.optsNum.optsSolv;
    else
        opts = struct('noIterations',40,'lambda',0.8,'Seppecher_red',1);
    end
           
    solveSquare   = [];

    if(~isfield(opts,'Seppecher_red'))
         Seppecher_red = 0;
    else
         Seppecher_red = opts.Seppecher_red;
    end
     
    Seppecher     = IsSeppecher(this);        
    %solveSquare = true;
    
    if(isfield(this.optsPhys,'theta'))
        theta = this.optsPhys.theta;
    elseif(~isempty(this.theta))
        theta = this.theta;
    else
        theta = pi/2;
    end
    
    if(Seppecher)            
        a      = [0;0];
        this.a = a;          
	end
      
    nParticles     = this.optsPhys.nParticles;            
    M              = this.IDC.M;  
    Ind            = this.IDC.Ind;
    IntSubArea     = this.IntSubArea;            
    
    IBB  = repmat(Ind.bound,2,1);      
    F    = false(M,1);   T       = true(M,1);        
    
    opts.optsNum    = this.optsNum;
    opts.optsPhys   = this.optsPhys;
    opts.configName = this.configName;    
     	            
	[res,~,Parameters] = DataStorage('ContactLineComputations',@SolveBinaryFluid,opts,[],[],{'optsNum_SubArea','optsNum_PlotArea','noIterations','Function'});                
        
    SetResults(res);

    if(IsSeppecher(this))
        Seppecher_red      = 2;
        opts.Seppecher_red = Seppecher_red;
        opts.lambda        = 0.6;    
        
        [res,~,Parameters] = DataStorage('ContactLineComputations',@SolveBinaryFluid,opts,[],[],{'optsNum_SubArea','optsNum_PlotArea','noIterations','Function'});    
        SetResults(res);
    end
	           
    this.filename               = Parameters.Filename;
    this.errors.errorIterations = res.errorHistory;        
    
    function SetResults(res)
        this.uv       = res.uv;
        this.mu       = res.mu;
        this.phi      = res.phi;
        this.p        = res.p;
        if(Seppecher_red == 1)
            this.a        = a;        
            this.deltaX   = res.deltaX;
            this.theta    = theta;    
        elseif(Seppecher_red == 2)
            this.a        = a;        
            this.deltaX   = res.deltaX;
            this.theta    = res.theta;    
        elseif(Seppecher)
            this.a        = res.a;        
            this.deltaX   = res.deltaX;
            this.theta    = theta;
            %this.theta    = res.theta;
        end
    end
        
    function [res,Results] = SolveBinaryFluid(opts,in)   
        
         if(Seppecher_red == 1)
             iGuess = GetInitialCondition(this,theta);
             iGuess = iGuess([3,5:end]);
         elseif(Seppecher_red == 2)     
              iGuess = [this.deltaX;this.theta;...
                        this.uv;this.phi;this.mu;this.p]; 
             %iGuess = GetInitialCondition(this,theta);
             %iGuess = iGuess(3:end);
         elseif(Seppecher)
             iGuess = GetInitialCondition(this,theta);
             iGuess = iGuess([1:3,5:end]);
         else
             iGuess = GetInitialCondition(this);        
         end
        %solveSquare = true;        

        solveSquare      = true;
        [vec,errHistory] = NewtonMethod(iGuess,@f,1e-6,opts.noIterations,opts.lambda);
        if(isempty(vec))            
            opts.optsNum.PhysArea
            opts.optsPhys
            opts.configName
            errHistory
            error('Newtonmethod failed to converge');                                    
        end

        solveSquare = false;
        optsFS = optimoptions(@fsolve,'Jacobian','on','Display','iter','Algorithm','levenberg-marquardt');
        vec = fsolve(@f,vec,optsFS);                    
        
        %[uv;phi;mu] 
        if(Seppecher_red == 1)
            res.deltaX = vec(1);            
            vec        = vec(2:end);
        elseif(Seppecher_red == 2)
            res.deltaX = vec(1);
            res.theta  = vec(2);
            vec        = vec(3:end);
        elseif(Seppecher)
            res.a      = vec(1:2);
            res.deltaX = vec(3);
            res.theta  = vec(4);
            vec        = vec(5:end);
        end
        res.uv     = vec([T;T;F;F;F]);
        res.phi    = vec([F;F;T;F;F]);
        res.mu     = vec([F;F;F;T;F]);
        res.p      = vec([F;F;F;F;T]);
        
        if(isfield(res,'deltaX'))
            Results.deltaX = res.deltaX;
        end
        if(isfield(res,'theta'))
            Results.theta = res.theta;
        end
        
        res.errorHistory = errHistory;
    end    
    function [v,A] = f(z)
                
        %[uv;phi;G;p]       
        if(Seppecher_red == 1)                        
            deltaX = z(1);            
            disp(['[deltaX,theta] = ',num2str(deltaX),' , ',num2str(theta*180/pi),'.']);        
            z      = z(2:end);            
        elseif(Seppecher_red == 2)
            deltaX = z(1);
            theta  = z(2);
            disp(['[deltaX,theta] = ',num2str(deltaX),' , ',num2str(theta*180/pi),'.']);        
            z      = z(3:end);                        
        elseif(Seppecher)
            a      = z(1:2);
            deltaX = z(3);
            %theta  = z(4);
            disp(['[a,deltaX,theta] = ',num2str(a(1)),' , ',num2str(a(2)),' , ',num2str(deltaX),' , ',num2str(theta*180/pi),'.']);
        
            z      = z(4:end);
        else
            a = []; deltaX = []; theta = [];                    
        end                

        uv  = z([T;T;F;F;F]);
        phi = z([F;F;T;F;F]);
        G   = z([F;F;F;T;F]);
        p   = z([F;F;F;F;T]);               
        
        %if(Seppecher_red)
%            a      = FindAB(this,phi,G,deltaX,a);
%        end

        [v_cont,A_cont] = Continuity(this,uv,phi,G,p);
        [v_mom,A_mom]   = Momentum(this,uv,phi,G,p);
        [v_G,A_G]       = PhasefieldEq(this,uv,phi,G);
        [v_mu,A_mu]     = ChemicalPotential(this,phi,G);
       
        A_particles              = zeros(1,5*M);
        A_particles([F;F;T;F;F]) = IntSubArea;
        v_particles              = IntSubArea*phi - nParticles;         

        if(Seppecher_red == 1)   
            A_cont      = [zeros(M,1),A_cont];
            A_mom       = [zeros(2*M,1),A_mom];
            A_G         = [zeros(M,1),A_G];
            A_mu        = [zeros(M,1),A_mu];
            A_particles = [0,A_particles];
            
            [v_SeppAdd,A_SeppAdd] = GetSeppecherConditions(this,uv,phi,G,p,deltaX,theta);
            v_SeppAdd = v_SeppAdd(3);
            A_SeppAdd = A_SeppAdd(3,[3,5:end]);
        elseif(Seppecher_red == 2)   
            A_cont      = [zeros(M,2),A_cont];
            A_mom       = [zeros(2*M,2),A_mom];
            A_G         = [zeros(M,2),A_G];
            A_mu        = [zeros(M,2),A_mu];
            A_particles = [0,0,A_particles];
            
            [v_SeppAdd,A_SeppAdd] = GetSeppecherConditions(this,uv,phi,G,p,deltaX,theta);
            v_SeppAdd   = v_SeppAdd(3:4);
            A_SeppAdd   = A_SeppAdd(3:4,[3:end]);            
        elseif(Seppecher)    
            A_cont      = [zeros(M,3),A_cont];
            A_mom       = [zeros(2*M,3),A_mom];
            A_G         = [zeros(M,3),A_G];
            A_mu        = [zeros(M,3),A_mu];
            A_particles = [0,0,0,A_particles];
            
            [v_SeppAdd,A_SeppAdd] = GetSeppecherConditions(this,uv,phi,G,p,deltaX,theta);
            v_SeppAdd = v_SeppAdd(1:3);
            A_SeppAdd = A_SeppAdd(1:3,[1:3,5:end]);
        end
        
        [v_mu(Ind.top|Ind.bottom),A_mu_BC] = GetPhiBC(this,phi,theta);          
        [v_mom(IBB),A_mom_BC]              = GetVelBC(this,uv,a,deltaX,theta);                    
        [v_G(Ind.bound),A_G_BC]            = GetChemPotBC(this,uv,phi,G,theta);
        
        if(Seppecher_red == 1)
            A_mu(Ind.top|Ind.bottom,:) = A_mu_BC(:,[3,5:end]);
            A_mom(IBB,:)               = A_mom_BC(:,[3,5:end]);
            A_G(Ind.bound,:)           = A_G_BC(:,[3,5:end]);
        elseif(Seppecher_red == 2)
            A_mu(Ind.top|Ind.bottom,:) = A_mu_BC(:,[3:end]);
            A_mom(IBB,:)               = A_mom_BC(:,[3:end]);
            A_G(Ind.bound,:)           = A_G_BC(:,[3:end]);
        elseif(Seppecher)
            A_mu(Ind.top|Ind.bottom,:) = A_mu_BC(:,[1:3,5:end]);
            A_mom(IBB,:)               = A_mom_BC(:,[1:3,5:end]);
            A_G(Ind.bound,:)           = A_G_BC(:,[1:3,5:end]);
        else
            A_mu(Ind.top|Ind.bottom,:) = A_mu_BC;
            A_mom(IBB,:)               = A_mom_BC;
            A_G(Ind.bound,:)           = A_G_BC;
        end
        
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
                
      %  A = A(2:end,2:end);
      %  v = v(2:end);
        
        %A = A([1,3:end],[1,3:end]);
        %v = v([1,3:end]);
    end    
    function DisplayError(error)
        
        if(Seppecher_red == 1)
            PrintErrorPos(error(1),'zero density at interface');        
            
            error = error(2:end);        
        elseif(Seppecher_red == 2)
            PrintErrorPos(error(1),'zero density at interface');        
            PrintErrorPos(error(2),'chemical potential at -inf');            
            error = error(3:end);                    
        elseif(Seppecher)
            PrintErrorPos(error(1),'consistent phase mass');
            PrintErrorPos(error(2),'consistent mass');
            PrintErrorPos(error(3),'zero density at interface');        
            %PrintErrorPos(error(4),'chemical potential at -inf');
        
            error = error(4:end);
        end
        PrintErrorPos(error([F;F;F;F;T]),'continuity equation',this.IDC.Pts);
        PrintErrorPos(error([T;F;F;F;F]),'y1-momentum equation',this.IDC.Pts);
        PrintErrorPos(error([F;T;F;F;F]),'y2-momentum equation',this.IDC.Pts);                            
        PrintErrorPos(error([F;F;F;T;F]),'G equation',this.IDC.Pts);
        PrintErrorPos(error([F;F;T;F;F]),'mu equation',this.IDC.Pts);
    end

end