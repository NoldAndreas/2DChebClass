function data = DDFT_DiffusionBox_2Phase(optsPhys,optsNum)
%************************************************************************* 
% data = DDFT_DiffusionBox_2Phase(optsPhys,optsNum)
%
% Set mu_s :=  = kBT*log(rho) + int(rho(r')*Phi2D(r-r'),dr') +...
%               + mu_HS(rho_s) + V_ext - mu
%
% Equilibrium:
% (EQ 1) mu_s     = 0  
% (EQ 2) int(rho) = Nparticles
% Dynamics: 
% (DYN 1) drho/dt = div(rho*grad(mu_s))
% (BC )    0      = n*grad(mu_s) 
% where n is the normal to the wall.
%*************************************************************************   
    disp(['** ',optsNum.DDFTCode,' **']);

    %************************************************
    %***************  Initialization ****************
    %************************************************
    close all;                    
    PhysArea    = optsNum.PhysArea;   
    N1 = PhysArea.N(1); N2 = PhysArea.N(2);              
    kBT          = optsPhys.kBT;     
    plotTimes    = optsNum.plotTimes;      
    if(isfield(optsPhys,'HSBulk'))
        HS_f       = str2func(optsPhys.HSBulk);  
    elseif(isfield(optsNum,'HSBulk'))
        HS_f       = str2func(optsNum.HSBulk); 
    else
        HS_f       = @ZeroMap;
    end                     
    nParticles = optsPhys.nParticlesS; 
    SubArea = optsNum.SubArea;                              
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    tic
	abox = Box(PhysArea);
    [Pts,Diff,Int,Ind,Interp] = abox.ComputeAll(optsNum.PlotArea); 
    opts = optsPhys; opts.optsNum = optsNum;
    convStruct     = DataStorage(['Box' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,abox);   
    Conv           = convStruct.Conv;
    
    Int_of_path                    = abox.IntFluxThroughDomain(100);
           
    subBox      = Box(SubArea);
    int         = subBox.IntFluxThroughDomain(100);
    IP          = abox.SubShapePts(subBox.Pts);
    Int_of_path = int*blkdiag(IP,IP);
            
    [Vext,Vext_grad]  = getVBackDVBack(Pts.y1_kv,Pts.y2_kv,optsPhys.V1);  
    
    %[Path,InterpPath,Int_of_path,Int_SubOnFull] = SubSpace(SubArea,...
%                @SpectralSpectral_Interpolation,Pts,Maps,'normal');        

    [~,rhoLiq_sat,mu_sat] = BulkSatValues(optsPhys,[0.01;0.6;-2]);
    TestAccuracyRectangle(abox,PhysArea,Conv,optsPhys,mu_sat,rhoLiq_sat,HS_f);%Pts,Maps,Interp
    
    I       = eye(N1*N2);     eyes    = [I I];
    
    t_preprocess = toc;
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    tic 
    x_ic    = fsolve(@f,zeros(1+N1*N2,1));        
    mu      = x_ic(1);
    x_ic    = x_ic(2:end);
    t_eqSol = toc;
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    tic
    mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([1;mM]));
    [outTimes,X_t] =  ode15s(@dx_dt,plotTimes,[0;x_ic],opts);        
    t_dynSol = toc;
    %************************************************
    %****************  Postprocess  ****************
    %************************************************        

    accFlux   = X_t(:,1);
    X_t       = X_t(:,2:end)';
    rho_t     = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
    
    flux_t    = zeros(2*N1*N2,length(plotTimes));
    for i = 1:length(plotTimes)
        flux_t(:,i) = GetFlux(X_t(:,i),plotTimes(i));
    end
    
%    Subspace = v2struct(Path,InterpPath,Int_of_path,Int_SubOnFull,accFlux);
    data     = v2struct(Conv,X_t,rho_t,...%Subspace
                        mu,flux_t,...
                        t_preprocess,t_eqSol,t_dynSol);      
    data.shape = abox;
    
    if(~isfield(optsNum,'savefileDDFT'))
        SaveToFile(optsNum.DDFTCode,v2struct(data,optsPhys,optsNum),getResultsPath());
    end                    
                    
	display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
    display(['Equilibrium Sol., Computation time (sec): ', num2str(t_eqSol)]);
    display(['Dynamics, Computation time (sec): ', num2str(t_dynSol)]);
    
    PlotDDFT(v2struct(optsPhys,optsNum,data));    
             
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function dxdt = dx_dt(t,x)
        
        x       = x(2:end);                

        mu_s     = GetExcessChemPotential(x,t,mu); 
        h_s      = Diff.grad*x - Vext_grad;
                
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        
        %Boundary Conditions: no flux at the walls        
        flux_dir        = Diff.grad*mu_s;
        dxdt(Ind.bound) = Ind.normal*flux_dir;           
             
        dxdt = [Int_of_path*GetFlux(x,t) ;dxdt];
    end

    function y = f(x)
        %solves for T*log*rho + Vext                
        mu_s         = x(1);
        x            = x(2:end);        
        rho_full     = exp((x-Vext)/kBT);
        y            = GetExcessChemPotential(x,0,mu_s); 
        y            = [Int*(rho_full) - nParticles;y];
    end

    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = GetExcessChemPotential(x,t,mu); 
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end

    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp((x-Vext)/kBT);
        mu_s  = x + Conv*rho_s - mu_offset + ...
                HS_f(rho_s,kBT) + ...
                getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys.V1);                
    end

end