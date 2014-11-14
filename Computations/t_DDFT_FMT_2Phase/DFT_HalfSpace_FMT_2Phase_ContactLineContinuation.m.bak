function data = DFT_HalfSpace_FMT_2Phase_ContactLineContinuation(optsPhys,optsNum)
%************************************************************************* 
% data = DDFT_DiffusionBox_2Phase_Sat(optsPhys,optsNum)
%
% Set mu_s :=  = kBT*log(rho) + int(rho(r')*Phi2D(r-r'),dr') +...
%               + mu_HS(rho_s) + V_ext - mu_sat
%
% Equilibrium:
% (EQ 1) mu_s = 0  
% Dynamics: 
% (DYN 1) drho/dt = div(rho*grad(mu_s))
% (BC )    0      = n*grad(mu_s) 
% where n is the normal to the wall.
%*************************************************************************            

    if(nargin == 0)
        Phys_Area = struct('N',[35,50],'L1',4,'L2',2,'y2wall',0.,...
                       'N2bound',16,'h',1,'L2_AD',2.);

        Plot_Area = struct('y1Min',-8,'y1Max',8,'N1',80,'N2',80,...
                           'y2Min',0.5,'y2Max',10);

        %Sub_Area  = Phys_Area;
        Sub_Area = struct('y1Min',-1,'y1Max',1,'N',[20,20],...
                          'y2Min',0,'y2Max',2);    

        ConvIn    = struct('L',[],'L2',1.,'N',[50,50]);

        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                           'N1disc',24,'N2disc',24);%'Ncircle',10,

        optsNum = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                         'plotTimes',0:0.1:6,...
                         'Conv',ConvIn,'FexNum',Fex_Num,...
                         'DDFTCode','DDFT_HalfSpace_FMT_2Phase_Sat');

        V1 = struct('V1DV1','Vext_Cart_9',...                          
                          'epsilon_w',0.8,'tau',1,...
                          'epsilon_w_end',0.0,...
                          'epsilon_w1',0,'w1_steepness',1);

        V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);	

        optsPhys = struct('V1',V1,'V2',V2,...                      
                          'kBT',0.7,...   
                          'Dmu',-0.05,'nSpecies',1,...
                          'sigmaS',1);
    end
    disp(['** ',optsNum.DDFTCode,' **']);
 
    %************************************************
    %***************  Initialization ****************
    %************************************************
    close all;
        
    PhysArea     = optsNum.PhysArea;   
    N1           = PhysArea.N(1); 
    N2           = PhysArea.N(2);              
    kBT          = optsPhys.kBT;    
    R            = optsPhys.sigmaS/2;
    plotTimes    = optsNum.plotTimes;      
    
	optsPhys.HSBulk = (['FexBulk_',optsNum.FexNum.Fex]);      
	getFex          = str2func(['Fex_',optsNum.FexNum.Fex]);
    nSpecies        = 1;
        
    SubArea        = optsNum.SubArea;
    Phi_r          = str2func(optsPhys.V2.V2DV2);
    Dmu            = optsPhys.Dmu;
    
	saveFigs = false;
       
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    tic        
    %[rhoGas_sat,rhoLiq_sat,mu_sat] = BulkPhaseDiagram(optsPhys);
    [rhoGas_sat,rhoLiq_sat,mu_sat] = BulkSatValues(optsPhys,[0.01;0.6;-2]);
    %GetCriticalPoint(optsPhys);

    optsPhys.mu_sat      = mu_sat;
    optsPhys.rhoGas_sat  = rhoGas_sat;
    optsPhys.rhoLiq_sat  = rhoLiq_sat;
    
	shape              = PhysArea;
    shape.Conv         = optsNum.Conv;
    HS                 = HalfSpace_FMT(shape,diag(optsPhys.sigmaS)/2);
    [Pts,Diff,Int,Ind] = HS.ComputeAll(optsNum.PlotArea);
    
    
    %****** Compute Convolution Matrix ***********
	opts.V2                            = optsPhys.V2;    
    opts.nSpecies                      = optsPhys.nSpecies;
    opts.optsNum.Conv                  = optsNum.Conv;    
	opts.optsNum.PhysArea              = optsNum.PhysArea;    
    [convStruct,recomputed,parameters] = DataStorage([class(HS) filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,HS);
    Conv                               = convStruct.Conv;
    
    if(recomputed)
        yPtsCheck          = [0 2 ; 0 PhysArea.y2wall ; -10 0 ; 20 10; 0 Pts.y2_kv(end-1); 0 inf];
        HS.TestConvolutionMatrix(yPtsCheck,@Phi);
    end
    
    %****** Compute FMT Matrices ***********
    tic
    fprintf(1,'Computing Fex matrices ...\n');   
        
    params.sigmaS  = optsPhys.sigmaS;
    params.FexNum  = optsNum.FexNum;
    params.PhysArea = optsNum.PhysArea;
    
    params.Pts     = HS.Pts;     
    params.Polar   = 'cart';      
    func           = str2func(['FexMatrices_',optsNum.FexNum.Fex]);                
    IntMatrFex     = DataStorage([class(HS) filesep func2str(func)],func,params,HS); %true 
    
    %********************************************
    
    Int_of_path = HS.IntFluxThroughDomain(100);
           
    subBox      = Box(SubArea);
    int         = subBox.IntFluxThroughDomain(100);
    IP          = HS.SubShapePts(subBox.Pts);
    Int_of_path = int*blkdiag(IP,IP);
            
    [Vext,Vext_grad]  = getVBackDVBack(Pts.y1_kv,Pts.y2_kv,optsPhys.V1);  
    figure('color','white');
    HS.doPlots(Vext+getVAdd(Pts.y1_kv,Pts.y2_kv,0,optsPhys.V1),'SC');
    I           = eye(N1*N2);  
    eyes        = [I I];
    t_preprocess = toc;
    
	%****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %****************************************************************        
    opts             = PhysArea;
    opts.optsPhys    = optsPhys;
    opts.FexNum      = optsNum.FexNum;
    opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,{'tau','epsilon_w_end'});     
    
    optss            = opts.optsPhys;
	optss.rho_iguess = rhoGas_sat;        
    optss.PhysArea   = Phys_Area;
                         
    %FMT_1DContinuation(HS,IntMatrFex,optss,opts.FexNum,Conv,'mu');    
    FMT_1DContinuation(HS,IntMatrFex,optss,opts.FexNum,Conv,'epw');    
        
    rho_ic1D         = Compute1DDensityProfile(opts,HS);  %rho_ic1D =  DataStorage('1DEquilibriumSolutions',@Compute1DDensityProfile,opts,[],true);     
    rho_ic           = repmat(rho_ic1D,HS.N1,1);    
    
    pause(2);
    close all;    
    
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    tic    
    fprintf('Solving for equilibrium condition...');                
    x_iguess  = kBT*log(rho_ic);    
    
    x_cont    = IterativeContinuation(@fEqCont_epw,40,0.2,[0;x_iguess],opts);
      
    n_cont     = size(x_cont,1);
    %x_ic      = DataStorage('EquilibriumSolutions',@ComputeEquilibriumCondition,opts,x_iguess);    
    par_cont   = x_cont(:,1);
    rho_cont   = exp((x_cont(:,2:end)'-Vext*ones(1,n_cont))/kBT);
    
    Gamma      = zeros(size(par_cont));
    markInf    = (Pts.y1_kv == inf);
    for i = 1:n_cont        
        Gamma(i) =  Int*(rho_cont(:,i)-repmat(rho_cont(markInf,i),N1,1));
    end
    
    figure('Color','white');
    for i = 1:n_cont
        rho = rho_cont(:,i);

        optDetails.nContours = [0.1,0.2,0.3,0.4,0.5,0.6];
        optDetails.clabel = true;
        
        subplot(1,2,1);        
        HS.doPlots(rho,'contour',optDetails);
        %HS.doPlots(rho,'SC');
        
        subplot(1,2,2);        
        plot(par_cont,Gamma,'linewidth',1.5); hold on;
        plot(par_cont(i),Gamma(i),'o','MarkerSize',10,'MarkerFace','g'); hold off;
        
        drawnow;
        pause(0.5);
    end
    
    t_eqSol   = toc;
    fprintf([num2str(t_eqSol),'s']);
%     %****************************************************************
%     %****************  Compute time-dependent solution   ************
%     %****************************************************************
%     tic
%     mM            = ones(N1*N2,1);
%     mM(Ind.bound) = 0;
%     opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([1;mM]));
%     [outTimes,X_t] =  ode15s(@dx_dt,plotTimes,[0;x_ic],opts);        
%     t_dynSol = toc;
    
    %************************************************
    %****************  Postprocess  ****************
    %************************************************        
    tic
    accFlux   = X_t(:,1);
    X_t       = X_t(:,2:end)';
    rho_t     = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
    
    flux_t    = zeros(2*N1*N2,length(plotTimes));
    for i = 1:length(plotTimes)
        flux_t(:,i) = GetFlux(X_t(:,i),plotTimes(i));
    end
    
    bulkData = v2struct(rhoGas_sat,rhoLiq_sat,mu_sat);
%    Subspace = v2struct(Path,InterpPath,Int_of_path,Int_SubOnFull,accFlux);
    data     = v2struct(X_t,rho_t,...%Subspace,...
                        flux_t,...
                        t_preprocess,t_eqSol,t_dynSol,bulkData);
	data.shape = HS;
    data.mu    = mu_sat + Dmu;
	
    
    if(~isfield(optsNum,'savefileDDFT'))
        SaveToFile(optsNum.DDFTCode,v2struct(data,optsPhys,optsNum),getResultsPath());
    end
    t_postprocess = toc;
    
    display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
    display(['Equilibrium, Computation time (sec): ', num2str(t_eqSol)]);
    display(['Dynamics, Computation time (sec): ', num2str(t_dynSol)]);
    display(['Postprocessor, Computation time (sec): ', num2str(t_postprocess)]);
    
    PlotDDFT(v2struct(optsPhys,optsNum,data));    
        
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************  
    
	function y = fEqCont_epw(x,epw)
        %solves for T*log*rho + Vext                           
        optsPhys.V1.epsilon_w1 = epw;        
        y                      = GetExcessChemPotential(x,0,mu_sat + Dmu);
    end    

    function y = fEqCont(x,mus)
        %solves for T*log*rho + Vext                        
        y            = GetExcessChemPotential(x,0,mu_sat + mus);                
    end
    
    function rho1D = Compute1DDensityProfile(opts,hs)   
        optss            = opts.optsPhys;
        optss.rho_iguess = rhoGas_sat;
        rho1D            = FMT_1D(HS,IntMatrFex,optss,opts.FexNum,Conv,true);
    end
    function x_ic = ComputeEquilibriumCondition(params,x_iguess)
        x_ic    = fsolve(@f,x_iguess);
    end

    function y = f(x)
        y  = GetExcessChemPotential(x,0,mu_sat + Dmu); 
    end

    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp(x/kBT);
        mu_s  = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end
        mu_s = mu_s + x + Conv*rho_s + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys.V1);        
    end   

    function dxdt = dx_dt(t,x)
        
        x       = x(2:end);     
        mu_s    = GetExcessChemPotential(x,t,mu_sat+Dmu); 
        
        h_s     = Diff.grad*x - Vext_grad;
                
        dxdt    = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        
        %Boundary Conditions: no flux at the walls        
        flux_dir        = Diff.grad*mu_s;
        dxdt(Ind.bound) = Ind.normal*flux_dir;          
        
        dxdt            = [Int_of_path*GetFlux(x,t) ;dxdt];
    end


    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = GetExcessChemPotential(x,t,mu_sat + Dmu); 
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end

    function z = Phi(y1,y2)
         z = Phi_r(sqrt(y1.^2 + y2.^2),optsPhys.V2);
    end       

end