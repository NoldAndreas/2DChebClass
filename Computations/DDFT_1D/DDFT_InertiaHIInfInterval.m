function [data,optsPhys,optsNum,optsPlot] = DDFT_InertiaHIInfInterval(optsPhys,optsNum,optsPlot)

    if(nargin == 0)
        [optsNum,optsPhys,optsPlot] = Default_DDFT_InertiaHIInfInterval();
    end
    
    AddPaths();
    
    close all;  
    disp(['** ',optsNum.DDFTCode,' **']);

    %************************************************
    %***************  Initialization ****************
    %************************************************    

    PhysArea    = optsNum.PhysArea;
    N = PhysArea.N;

    
    kBT         = optsPhys.kBT; 
    nParticlesS = optsPhys.nParticlesS;            
    nSpecies=length(nParticlesS);
    
    mS = optsPhys.mS;
    gammaS = optsPhys.gammaS;
    
    gamma  = bsxfun(@times,gammaS',ones(N,nSpecies));
    m  = bsxfun(@times,mS',ones(N,nSpecies));

    mInv  = m.^(-1);
    
    plotTimes   = optsNum.plotTimes;
    
    Fex_Num = optsNum.FexNum;
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************ 
    
    aLine = InfSpectralLine(PhysArea);

    [Pts,Diff,Int,Ind,~] = aLine.ComputeAll(optsNum.PlotArea); 

    yS               = repmat(Pts.y,1,nSpecies);
    
    
    if(strcmp(Fex_Num.Fex,'Percus'))
        optsFMT        = Fex_Num;
        optsFMT.sigma  = optsPhys.sigmaS;
        IntMatrFex     = aLine.ComputeFMTMatrices(optsFMT);
        getFex         = str2func([Fex_Num.Fex]);
        getConv        = @zeroFunction;
        convStruct     = [];
    elseif(strcmp(Fex_Num.Fex,'Meanfield'))
        getFex         = @zeroFunction;
        IntMatrFex     = [];
        opts           = optsPhys; 
        opts.optsNum   = optsNum;
        opts.FexNum    = Fex_Num;
        convStruct     = DataStorage(['Interval' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,aLine,false);
        getConv        = @Fex_Meanfield;
    elseif(strcmp(Fex_Num.Fex,'Perturbation'))
        optsFMT        = Fex_Num;
        optsFMT.sigma  = optsPhys.sigmaS;
        IntMatrFex     = aLine.ComputeFMTMatrices(optsFMT);
        getFex         = str2func([Fex_Num.Fex]);
        opts           = optsPhys; 
        opts.optsNum   = optsNum;
        opts.FexNum    = Fex_Num;
        convStruct     = DataStorage(['Interval' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,aLine,false);
        getConv        = @Fex_Meanfield;
    end
    
    opts.optsPhys = optsPhys;
    opts.optsNum = optsNum;
    HIStruct = DataStorage(['Interval' filesep 'HIMatrices'],@HIMatrices,opts,aLine,false);
       
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    [Vext,Vext_grad]  = getVBackDVBack1D(yS,0,optsPhys.V1);

    % add in VAdd here
    x0 = getInitialGuess(0);    
    mu0 = zeros(1,nSpecies);
    x0mu0 = [mu0;x0];
        
    x_ic    = fsolve(@f,x0mu0);
    mu      = x_ic(1,:);
    mu      = repmat(mu,N,1);
    x_ic    = x_ic(2:end,:);
    v_ic    = zeros(size(x_ic));

    xv_ic   = [x_ic ; v_ic];
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    mMx            = ones(N,1);
    mMv            = ones(N,1);
    %mMv(Ind.bound) = 0;
    mM             = [mMx;mMv];
    mM = repmat(mM,nSpecies,1);
    opts = odeset('RelTol',10^-6,'AbsTol',10^-6,'Mass',diag(mM));
    [~,XV_t] =  ode15s(@dxv_dt,plotTimes,xv_ic,opts);        

    XV_t = XV_t.';
    XV_t = reshape(XV_t,[],nSpecies,size(XV_t,2));
    
    X_t  = XV_t(1:end/2,:,:);
    
    rho_t     = zeros(N,nSpecies,length(plotTimes));
    flux_t    = XV_t(end/2+1:end,:,:);
    V_t       = zeros(N,nSpecies,length(plotTimes));
    for i = 1:length(plotTimes)
        rho_t(:,:,i)  = exp((X_t(:,:,i)-Vext)/kBT);
        V_t(:,:,i)    = Vext + getVAddDVAdd1D(yS,plotTimes(i),optsPhys.V1);
    end
    
    data       = v2struct(IntMatrFex,convStruct,HIStruct,X_t,rho_t,mu,flux_t,V_t);
    data.shape = aLine;

    if(~isfield(optsNum,'doPlots') ...
            || (isfield(optsNum,'doPlots') && optsNum.doPlots) )
        PlotDDFT1D(v2struct(optsPhys,optsNum,data));       
    end
    
    
   %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function dxvdt = dxv_dt(t,x)

        xv = reshape(x,[],nSpecies);
       
        x = xv(1:end/2,:);
        v = xv(end/2+1:end,:);
                
        mu_s = GetChemicalPotential(x,t,mu);
        
        h_s   = Diff.Dy*x - Vext_grad;
        
        dxdt  = - kBT*Diff.Dy*v - h_s.*v;  
        
        dvdt  = - gamma.*v - v.*(Diff.Dy*v) - mInv.*(Diff.Dy*mu_s);
        
        rho_s = exp((x-Vext)/kBT);
        HI = ComputeHI(rho_s,v,HIStruct);
        
        dvdt = dvdt - gamma.*HI;
        
        dxdt(Ind.bound,:) = 0;
        dvdt(Ind.bound,:) = 0;

        dxvdt = [dxdt;dvdt];
        
        dxvdt = dxvdt(:);         
 
    end
    
    function y = f(x)
        %solves for T*log*rho + Vext                
        mu_s         = x(1,:);
        mu_s         = repmat(mu_s,N,1);
        x            = x(2:end,:);        
        rho_full     = exp((x-Vext)/kBT);
        y            = GetChemicalPotential(x,0,mu_s);
        y            = [Int*rho_full - nParticlesS.';y];
    end
    
    function mu_s = GetChemicalPotential(x,t,mu)
        rho_s = exp((x-Vext)/kBT); 
        mu_s  = getFex(rho_s,IntMatrFex,kBT);
        VAdd  = getVAddDVAdd1D(yS,t,optsPhys.V1);
        ConvRho_s = getConv(rho_s,convStruct,kBT);
        mu_s  = mu_s + x + ConvRho_s - mu + VAdd;
    end

    function x0 = getInitialGuess(VAdd)
        rhoInit = exp(-(Vext+VAdd)/kBT);
        normalization = repmat( Int*rhoInit./nParticlesS' , size(rhoInit,1) ,1);
        x0=-kBT*log(normalization) - VAdd;
    end

    function mu = zeroFunction(rho_s,~,~)
        mu = zeros(size(rho_s));
    end

end