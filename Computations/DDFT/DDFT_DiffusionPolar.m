function output = DDFT_DiffusionPolar(optsPhys,optsNum)
%************************************************
%output = DDFT_DiffusionPolar(optsPhys,optsNum)
%At Equilibrium:
% (1 EQUIL) 0 = kBT*log(rho) + int(rho*Phi2D) - mu + Vext
% (2 EQUIL) n = int(rho)
%Dynamics:
% define    mu_s = kBT*log(rho) + int(rho*Phi2D) - mu + Vext
% (1 DYN) drho/dt = div(rho*grad(mu_s))
% (BC1)   0       = n*(grad(mu_s))
% if full domain, 
% (BC Inf) rho_ic = rho(inf)
%************************************************

    if(nargin == 0)
         [optsNum,optsPhys] = Default1_DDFT_DiffusionPolarDisk();
    end
    
    disp(optsNum.DDFTCode);
    
    %************************************************
    %****************  Default values****************
    %************************************************    
    close all;
    kBT = optsPhys.kBT; nParticles = optsPhys.nParticlesS;        
    [N1,N2,PhysArea,h1,plotTimes] = LoadNumData(optsNum);
    
    %************************************************
    %***************  Preprocessing  ****************
    %************************************************    
    tic
    R = PhysArea.y1Max;
    N = [N1;N2];
    DC = Disc(v2struct(R,N));
    [Pts,Diff,Int,Ind,Interp] = DC.ComputeAll(optsNum.PlotArea);	
    
    opts = optsPhys;    opts.optsNum = optsNum;
	convStruct        = DataStorage(['Disc' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,DC);    
    Conv              = convStruct.Conv;
    toc
    
    [Vext,Vext_grad]  = getVBackDVBack(Pts.y1_kv,Pts.y2_kv,optsPhys.V1);
    
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    x_ic    = fsolve(@f,zeros(1+N1*N2,1));    
    mu      = x_ic(1);
    x_ic    = x_ic(2:end);    
    
    DC.doPlots(exp((x_ic-Vext)/kBT));  
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    mM           = ones(N1*N2,1);
    mM(Ind.outR) = 0;    
    opts   = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
    [outTimes,X_t] = ode15s(@dx_dt,plotTimes,x_ic,opts);
    
    %************************************************
    %****************  Postprocess  *****************
    %************************************************
    X_t    = X_t';
    rho_t  = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
 
    flux_t    = zeros(2*N1*N2,length(plotTimes));
    for i = 1:length(plotTimes)
        flux_t(:,i) = GetFlux(X_t(:,i),plotTimes(i));
    end    
    
    fourier    = 'fourier';   
    data       = v2struct(Conv,X_t,rho_t,mu,flux_t,fourier);    
    data.shape = DC;

    if(~isfield(optsNum,'savefileDDFT'))
        SaveToFile(optsNum.DDFTCode,v2struct(data,optsPhys,optsNum),getResultsPath());
    end            
        
    PlotDDFT(v2struct(optsPhys,optsNum,data));   

    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function dxdt = dx_dt(t,x)
        
        rho_s       = exp((x-Vext)/kBT);
        x_full =x;
        
        I = eye(N1*N2);

        mu_s  = GetChemicalPotential(x,t); %mu_s  = x + Conv*rho_s - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys);
        mu_s(Pts.y1_kv==inf) = 0;
        
        h_s                                  = Diff.grad*x_full - Vext_grad;
        h_s([Pts.y1_kv==inf;false(N1*N2,1)]) = 0; %here, we have assumed that grad(mu) converges fast enough
        
        dxdt             = kBT*Diff.Lap*mu_s + [I I]*(h_s.*(Diff.grad*mu_s));  
        
        %Boundary Conditions at infinity
        if(PhysArea.y1Max == inf)
            dxdt(Ind.outR)  = x(Ind.outR) - x_ic(Ind.outR);   
        else
            flux_dir        = Diff.grad*mu_s;
            dxdt(Ind.outR)  = Ind.normalOutR*flux_dir;             
        end
        
    end    
    
    function y = f(x)
        %solves for T*log*rho + Vext      
        mu           = x(1);
        x            = x(2:end);                
        rho_full     = exp((x-Vext)/kBT);
        y            = GetChemicalPotential(x,0); %y = x + Conv*rho_full - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,0,optsPhys);        
        y            = [Int*rho_full - nParticles;y]; 
    end

    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s = GetChemicalPotential(x,t); %mu_s  = x + Conv*rho_s - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys);   
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end

    function mu_s = GetChemicalPotential(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = x + Conv*rho_s - mu + getVAdd(Pts.y1_kv,Pts.y2_kv,t,optsPhys.V1);
    end
    
end

