function data = DDFT_DiffusionHalfSpace(optsPhys,optsNum)
%************************************************************************
%  define mu_s = kBT*log(rho) + int(rho(r')*Phi2D(r-r'),dr') + V_ext - mu
% Equilibrium:
% (EQ 1)  mu_s      = 0
% (EQ 2)  int(rho)  = noParticles
% Dynamics:
% (DYN 1) drho/dt = div(rho*grad(mu_s))
% (BC 1)  n*grad(rho) = 0
%************************************************************************
    if(nargin == 0)
        [optsNum,optsPhys] = Default_HalfSpace_DDFT_DiffusionClosedBox();
    end

    disp(optsNum.DDFTCode); 
    %************************************************
    %***************  Initialization ****************
    %************************************************
    close all;
    
    kBT         = optsPhys.kBT;
    mu          = optsPhys.mu;    
    PhysArea    = optsNum.PhysArea;
    PlotArea    = optsNum.PlotArea;
    optsNum.plotTimes = (0:(optsNum.TN-1))*optsNum.Tmax/(optsNum.TN-1);
    N1          = PhysArea.N(1);  N2 = PhysArea.N(2);    
    Phi_r       = str2func(optsPhys.V2.V2DV2);

    %************************************************
    %****************  Preprocess  ****************
    %************************************************  
    shape        = PhysArea;
    shape.Conv   = optsNum.Conv;
   
    HS                  = HalfSpace(shape);
    [Pts,Diff,Int,Ind]  = HS.ComputeAll(PlotArea);    
    
	opts              = optsPhys;    
    opts.optsNum      = optsNum;
	convStruct        = DataStorage(['HalfSpace' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,HS);    
    Conv              = convStruct.Conv;
    %Conv             = HS.ComputeConvolutionMatrix(@Phi,true,yPtsCheck);        
    yPtsCheck         = [20 291; 0 2 ; 0 PhysArea.y2Min ; -10 0];
    HS.TestConvolutionMatrix(yPtsCheck,@Phi);
                                
    I       = eye(N1*N2);            
    eyes    = [I I];
    
    [Vext,Vext_grad]  = getVBackDVBack(Pts.y1_kv,Pts.y2_kv,optsPhys.V1);    
    
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************        
    x_ic    = fsolve(@f,zeros(N1*N2,1));
    rho_ic  = exp((x_ic-Vext)/kBT);
    HS.doPlots(rho_ic,'SC');
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    [outTimes,X_t] =  ode15s(@dx_dt,optsNum.plotTimes,x_ic,opts);        

    %************************************************
    %****************  Postprocess  ****************
    %************************************************        
    
    X_t       = X_t';
    rho_t     = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
    
    flux_t    = zeros(2*N1*N2,length(optsNum.plotTimes));
    for i = 1:length(optsNum.plotTimes)
        flux_t(:,i) = GetFlux(X_t(:,i),optsNum.plotTimes(i));
    end
        
    data       = v2struct(X_t,rho_t,mu,flux_t);
    data.shape = HS;
    
    if(~isfield(optsNum,'savefileDDFT'))
        SaveToFile(optsNum.DDFTCode,v2struct(data,optsPhys,optsNum),getResultsPath());
    end            
    
    PlotDDFT(v2struct(optsPhys,optsNum,data));                 
        
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
	function y = f(x)
        %solves for T*log*rho + Vext                        
        rho_full     = exp((x-Vext)/kBT);
        y            = x + Conv*rho_full - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,0,optsPhys.V1);                        
    end

    function dxdt = dx_dt(t,x)
        rho_s   = exp((x-Vext)/kBT);                
        
        mu_s  = x + Conv*rho_s - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys.V1);
        mu_s(Pts.y1_kv == inf) = 0;                
        mu_s(Pts.y1_kv == -inf) = 0;                
        
        h_s   = Diff.grad*x - Vext_grad;
        h_s([Pts.y1_kv==inf | Pts.y1_kv==-inf;false(N1*N2,1)]) = 0; %here, we have assumed that grad(mu) converges fast enough                                      
                
        dxdt            = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        
        %Boundary Conditions: no flux at the walls        
        flux_dir         = Diff.grad*mu_s;
        dxdt(Ind.bottom) = Ind.normalBottom*flux_dir;           
        dxdt(Ind.top)    = x(Ind.top) - x_ic(Ind.top);        
        dxdt(Ind.right)  = x(Ind.right) - x_ic(Ind.right);
        dxdt(Ind.left)   = x(Ind.left) - x_ic(Ind.left);                        
    end



    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = x + Conv*rho_s - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys.V1);                
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end

    function z = Phi(y1,y2)
         z = Phi_r(sqrt(y1.^2 + y2.^2),optsPhys.V2);
    end            


end