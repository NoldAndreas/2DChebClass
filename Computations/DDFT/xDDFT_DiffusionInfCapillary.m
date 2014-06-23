function data = DDFT_DiffusionInfCapillary(optsPhys,optsNum)
%************************************************
%  define mu_s = kBT*log(rho) + int(rho(r')*Phi2D(r-r'),dr') + V_ext - mu
% Equilibrium:
% (EQ 1)  mu_s      = 0
% (EQ 2)  int(rho)  = noParticles
% Dynamics:
% (DYN 1) drho/dt = div(rho*grad(mu_s))
% (BC 1)  n*grad(rho) = 0
%************************************************
    if(nargin == 0)
        [optsNum,optsPhys] = Default_InfiniteCapillary_DDFT_DiffusionClosedBox();
    end

    disp(optsNum.DDFTCode);
 
    %************************************************
    %***************  Initialization ****************
    %************************************************
    close all;
    
    kBT         = optsPhys.kBT;
    nParticles  = optsPhys.nParticlesS;    
    PhysArea    = optsNum.PhysArea;
    PlotArea    = optsNum.PlotArea;
    Phi_r       = str2func(optsPhys.V2.V2DV2);
    optsNum.plotTimes = (0:(optsNum.TN-1))*optsNum.Tmax/(optsNum.TN-1);
    N1          = PhysArea.N(1);  N2 = PhysArea.N(2);    

    %************************************************
    %****************  Preprocess  ****************
    %************************************************  
    shape        = PhysArea;
    shape.Conv   = optsNum.Conv;
   
    IC                         = InfCapillary(shape);
    [Pts,Diff,Int,Ind,Interp]  = IC.ComputeAll(PlotArea);
    
    yPtsCheck                  = [0 0 ; 0 PhysArea.y2Min ; -10 0 ; 20 0];
    
	opts = optsPhys;    opts.optsNum = optsNum;
	convStruct        = DataStorage(['InfCapillary' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,IC);    
    Conv              = convStruct.Conv;    
    IC.TestConvolutionMatrix(yPtsCheck,@Phi);
   
    subArea        = Box(optsNum.SubArea);
    IP             = IC.SubShapePts(subArea.Pts);
    %Int_SubOnFull = subArea.ComputeIntegrationVector()*IP;
        
    Int_of_path   = subArea.IntFluxThroughDomain(100)*blkdiag(IP,IP);
                                
    I       = eye(N1*N2);            
    eyes    = [I I];
    
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    [Vext,Vext_grad]  = getVBackDVBack(Pts.y1_kv,Pts.y2_kv,optsPhys.V1);
    x_ic    = fsolve(@f,zeros(1+N1*N2,1));
    mu      = x_ic(1);
    x_ic    = x_ic(2:end);
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([1;mM]));
    [outTimes,X_t] =  ode15s(@dx_dt,optsNum.plotTimes,[0;x_ic],opts);        

    %************************************************
    %****************  Postprocess  ****************
    %************************************************        

    accFlux   = X_t(:,1);
    X_t       = X_t(:,2:end)';
    rho_t     = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
    
    flux_t    = zeros(2*N1*N2,length(optsNum.plotTimes));
    for i = 1:length(optsNum.plotTimes)
        flux_t(:,i) = GetFlux(X_t(:,i),optsNum.plotTimes(i));
    end
    
    Subspace = v2struct(subArea,accFlux);%Path,InterpPath,Int_of_path,Int_SubOnFull
    data     = v2struct(X_t,rho_t,Subspace,mu,flux_t);
    data.shape = IC;
    
    if(~isfield(optsNum,'savefileDDFT'))
        SaveToFile(optsNum.DDFTCode,v2struct(data,optsPhys,optsNum),getResultsPath());
    end            
    
    PlotDDFT(v2struct(optsPhys,optsNum,data));                 
        
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function dxdt = dx_dt(t,x)
        
        x       = x(2:end);        
        rho_s   = exp((x-Vext)/kBT);                
        
        mu_s  = x + Conv*rho_s - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys.V1);
        mu_s(Pts.y1_kv == inf) = 0;                
        mu_s(Pts.y1_kv == -inf) = 0;                
        
        h_s   = Diff.grad*x - Vext_grad;
        h_s([Pts.y1_kv==inf | Pts.y1_kv==-inf;false(N1*N2,1)]) = 0; %here, we have assumed that grad(mu) converges fast enough                                      
                
        dxdt            = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        
        %Boundary Conditions: no flux at the walls        
        flux_dir        = Diff.grad*mu_s;
        dxdt(Ind.bound) = Ind.normal*flux_dir;           
        %dxdt(Ind.bound) = Ind.normal*flux_dir;     
        dxdt(Ind.right) = x(Ind.right) - x_ic(Ind.right);
        dxdt(Ind.left)  = x(Ind.left) - x_ic(Ind.left);
                
        dxdt = [Int_of_path*GetFlux(x,t) ;dxdt];
    end
    function y = f(x)
        %solves for T*log*rho + Vext                
        mu_s         = x(1);
        x            = x(2:end);        
        rho_full     = exp((x-Vext)/kBT);
        y            = x + Conv*rho_full - mu_s + getVAdd(Pts.y1_kv,Pts.y2_kv,0,optsPhys.V1);                
        y            = [Int*rho_full - nParticles;y];
    end
    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = x + Conv*rho_s - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys.V1);                
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end
    function z = Phi(r)%y1,y2)
         z = Phi_r(r,optsPhys.V2);%sqrt(y1.^2 + y2.^2)
    end            


end