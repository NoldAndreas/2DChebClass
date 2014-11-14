function data = FMT_2D_HardWall(optsPhys,optsNum,optsPlot)
%************************************************************************* 
% define
%   mu_s_i := kBT*log(rho_i) + sum( int(rho_j(r')*Phi2D(r-r'),dr') , 
%               j = 1..N) +mu_HS(rho_i) + V_ext_i - mu_i
% Equilibrium, (solve for each species i)
%  (EQ i,1) mu_s_i     = 0
%  (EQ i,2) int(rho_i) = NoParticles_i
% Dynamics: 
%*************************************************************************   
    if(nargin == 0)
        [optsNum,optsPhys,optsPlot] = Test2FMT_DDFT_DiffusionHalfSpace_HardWall();       
    end    
    saveFigs = false;

    close all;  
    disp(['** ',optsNum.DDFTCode,' **']);

    %************************************************
    %***************  Initialization ****************
    %************************************************
    PhysArea  = optsNum.PhysArea;           
    R         = diag(optsPhys.sigmaS)/2;    
    kBT       = optsPhys.kBT;      
    eta       = optsPhys.eta;
    nSpecies  = 1;
    fBulk     = str2func(['FexBulk_',optsNum.FexNum.Fex]);  
    getFex    = str2func(['Fex_',optsNum.FexNum.Fex]);
    
	rhoBulk   = eta*6/pi;
    mu        = kBT*log(rhoBulk) + fBulk(rhoBulk,kBT);
        
    N1        = PhysArea.N(1);   N2 = PhysArea.N(2);
    plotTimes = 0:(optsNum.Tmax/optsNum.TN):optsNum.Tmax;
	I         = eye(N1*N2);   eyes    = repmat(I,1,2);
        
    HS        = HalfSpace_FMT(PhysArea,R,optsNum.Accuracy_Averaging);    
    [Pts,Diff,Int,Ind,Interp] = HS.ComputeAll(optsNum.PlotArea);    
    
    y1S     = repmat(Pts.y1_kv,1,nSpecies); 
    y2S     = repmat(Pts.y2_kv,1,nSpecies);
    [Vext,Vext_grad]  = getVBackDVBack(y1S,y2S,optsPhys.V1);       

    %************************************************
    %****************  Preprocess  ******************
    %************************************************
    
    tic
    fprintf(1,'Computing Fex matrices ...\n');   
    
    params         = optsPhys;
    params         = rmfield(params,'V1');
    params.FexNum  = optsNum.FexNum;
    params.PhysArea = optsNum.PhysArea;
    
    params.Pts     = HS.Pts;     
    params.Polar   = 'cart';      
    func           = str2func(['FexMatrices_',optsNum.FexNum.Fex]);            
    params         = rmfield(params,'eta');
    
    IntMatrFex_2D  = DataStorage('FexData',func,params,HS); %true 
    
    
    %****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %****************************************************************        

    rho_ic1D       = FMT_1D_HardWall(HS,IntMatrFex_2D,optsPhys,optsNum);
    rho_ic         = repmat(rho_ic1D,HS.N1,1);
    x_ic           = kBT*log(rho_ic);
    
    pause(2);
    close all;       
    
    %****************************************************************
    %**************** Solve for equilibrium 2D condition   **********
    %****************************************************************        
    
    disp('Compute 2D profile..');        
        
   %x_ic           = zeros(size(x_ic));    
    icParams       = v2struct(optsPhys,optsNum);
    %icParams.optsPhys.V1 = rmfield(icParams.optsPhys.V1,'grav_end');
    %icParams.optsNum     = rmfield(icParams.optsNum,'Tmax');
    %icParams.optsNum     = rmfield(icParams.optsNum,'TN');
    
    icParams.fsolveOpts = optimset('MaxFunEvals',2000000,'MaxIter',200000);    
    x_ic           = DataStorage('EquilibriumData',@ComputeEquilibrium,icParams,x_ic);         	
    rho_ic         = exp(x_ic/kBT);
        
    figure
    
	HS.plot(rho_ic,'SC');
	%PlotRosenfeldFMT_AverageDensities(HS,IntMatrFex(1),rho_ic);        
    
    %****************************************************************
    %**************** Solve for dynamic 3D condition   **********
    %****************************************************************        
    disp('Compute dynamic 2D evolution ..');        
    mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
    [outTimes,X_t] =  ode15s(@dx_dt,[0 max(plotTimes)],x_ic,opts);                
        
    plotTimes = outTimes;
    X_t       = X_t';
    rho_t     = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
    
    flux_t    = zeros(2*N1*N2,length(plotTimes));
    for i = 1:length(plotTimes)
        flux_t(:,i) = GetFlux(X_t(:,i),plotTimes(i));
    end
    
    data     = v2struct(X_t,rho_t,mu,flux_t);%Pts,Diff,Int,Ind,Interp
    data.shape = HS;
    optsNum.plotTimes = plotTimes; 
    
    plotData = v2struct(optsPhys,optsNum,data);
    
    SaveToFile(optsNum.DDFTCode,plotData,getResultsPath());
    PlotDDFT(plotData); 
    
    
        
    %***********************************************************
    %**************** Functions ********************************
    %***********************************************************
    
    function y = ComputeEquilibrium(params,y0)
        fprintf(1,'Computing initial condition ...');        
        y          = fsolve(@f,y0,params.fsolveOpts);     
    end
    function y = f(x)
        %solves for T*log*rho + Vext                        
        y            = GetExcessChemPotential(x,0,mu);         
        y            = y(:);
    end

    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp(x/kBT);                
        mu_s  = getFex(rho_s,IntMatrFex_2D,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end        
        mu_s = mu_s + x + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys.V1);  %HS_f(rho_s,kBT)
    end

    function dxdt = dx_dt(t,x)
       % x       = reshape(x,N1*N2,nSpecies);
        
        mu_s     = GetExcessChemPotential(x,t,mu);
        mu_s((Pts.y1_kv==inf) | (Pts.y1_kv==-inf) | (Pts.y2_kv==inf),:) = 0;
        h_s      = Diff.grad*x - Vext_grad;
        h_s([Pts.y2_kv==inf;Pts.y2_kv==inf],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==inf;false(N1*N2,1)],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==-inf;false(N1*N2,1)],:) = 0; %here, we have assumed that grad(mu) converges fast enough
                
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
                
        %Boundary Conditions at infinity                
        flux_dir           = Diff.grad*mu_s;                
        dxdt(Ind.bottom,:) = Ind.normalBottom*flux_dir;             
               
        dxdt(Ind.top,:)    = x(Ind.top,:) - x_ic(Ind.top,:);   
        dxdt(Ind.right,:)  = x(Ind.right,:) - x_ic(Ind.right,:);   
        dxdt(Ind.left,:)   = x(Ind.left,:) - x_ic(Ind.left,:);   
        dxdt = dxdt(:);      
    end

    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = GetExcessChemPotential(x,t,mu); 
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end


end