function [data,optsPhys,optsNum,optsPlot] = DDFT_DiffusionSphericalInfInterval(optsPhys,optsNum,optsPlot)

    if(nargin == 0)
        %[optsNum,optsPhys,optsPlot] = Default_DDFT_DiffusionSphericalInfInterval();
        [optsNum,optsPhys,optsPlot] = Default_DDFT_DiffusionHSSphericalInfInterval();
       
    end
    
    AddPaths();
    
    close all;  
    disp(['** ',optsNum.DDFTCode,' **']);

    %************************************************
    %***************  Initialization ****************
    %************************************************    

    PhysArea    = optsNum.PhysArea;
    N = PhysArea.N;

    if(mod(N,2)==1)
        N=N+1;
    end
    
    kBT         = optsPhys.kBT; 
    nParticlesS = optsPhys.nParticlesS;            
    nSpecies=length(nParticlesS);
    
    mS = optsPhys.mS;
    gammaS = optsPhys.gammaS;
    gamma  = bsxfun(@times,gammaS',ones(N,nSpecies));
    m  = bsxfun(@times,mS',ones(N,nSpecies));
    D0 = 1./(gamma.*m);
    
    plotTimes   = optsNum.plotTimes;
    
    Fex_Num = optsNum.FexNum;
            
    %************************************************
    %****************  Preprocess  ****************
    %************************************************ 
    
    aLine = InfSpectralLineSpherical(PhysArea);
    
    [Pts,Diff,~,Ind,~] = aLine.ComputeAll(optsNum.PlotArea); 
    
    % SORT THE HARD CODING HERE
    geom.yMin = 0;
    geom.yMax = 20;
    IntNP = aLine.IntegrateRegion(@NPweight,geom);
       
    yS               = repmat(Pts.y,1,nSpecies);
    
    if(strcmp(Fex_Num.Fex,'FMT'))
        opts           = optsPhys;
        IntMatrFex     = DataStorage(['Interval' filesep 'FexMatrices_FMTRosenfeldSpherical'],@FexMatrices_FMTRosenfeldSpherical,opts,aLine,false);
        getFex         = @Fex_FMTRosenfeldSpherical;
        getConv        = @zeroFunction;
        convStruct     = [];
    elseif(strcmp(Fex_Num.Fex,'Meanfield'))
        getFex         = @zeroFunction;
        IntMatrFex     = [];
        opts           = optsPhys; 
        opts.optsNum   = optsNum;
        opts.FexNum    = Fex_Num;
        convStruct     = DataStorage(['Interval' filesep 'FexMatrices_Meanfield_Spherical'],@FexMatrices_Meanfield_Spherical,opts,aLine,false);
        getConv        = @Fex_Meanfield;
    elseif(strcmp(Fex_Num.Fex,'Perturbation'))
        opts           = optsPhys;
        IntMatrFex     = DataStorage(['Interval' filesep 'FexMatrices_FMTRosenfeldSpherical'],@FexMatrices_FMTRosenfeldSpherical,opts,aLine,false);
        getFex         = @Fex_FMTRosenfeldSpherical;
        opts           = optsPhys; 
        opts.optsNum   = optsNum;
        opts.FexNum    = Fex_Num;
        convStruct     = DataStorage(['Interval' filesep 'FexMatrices_Meanfield_Spherical'],@FexMatrices_Meanfield_Spherical,opts,aLine,false);
        getConv        = @Fex_Meanfield;
    end
          
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    [Vext,Vext_grad]  = getVBackDVBack1D(yS,0,optsPhys.V1);
    
    % add in VAdd here
    x0 = getInitialGuess(0);   
    x0 = cut(x0);
    
    mu0 = zeros(1,nSpecies);
    x0mu0 = [mu0;x0];
        
    x_ic    = fsolve(@f,x0mu0);
    mu      = x_ic(1,:);
    mu      = repmat(mu,N,1);
    x_ic    = x_ic(2:end,:);
   
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    mM = ones(N/2,1);
    mM = repmat(mM,nSpecies,1);
    opts = odeset('RelTol',10^-6,'AbsTol',10^-6,'Mass',diag(mM));
    [~,X_t] =  ode15s(@dx_dt,plotTimes,x_ic,opts);        

    X_t = X_t.';
    X_t = reshape(X_t,[],nSpecies,size(X_t,2));
    
    rho_t     = zeros(N,nSpecies,length(plotTimes));
    flux_t    = zeros(N,nSpecies,length(plotTimes));
    V_t       = zeros(N,nSpecies,length(plotTimes));
    for i = 1:length(plotTimes)
        x_t = mirror(X_t(:,:,i));
        rho_t(:,:,i)  = exp((x_t-Vext)/kBT);
        flux_t(:,:,i) = GetFlux(x_t,plotTimes(i));
        V_t(:,:,i)    = Vext + getVAddDVAdd1D(yS,plotTimes(i),optsPhys.V1);
    end
    
    data       = v2struct(IntMatrFex,convStruct,X_t,rho_t,mu,flux_t,V_t);
    data.shape = aLine;
    data.shape.Int = IntNP;
    if(~isfield(optsNum,'doPlots') ...
            || (isfield(optsNum,'doPlots') && optsNum.doPlots) )
        PlotDDFT1D(v2struct(optsPhys,optsNum,data));       
    end
    
   %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function dxdt = dx_dt(t,x)
        
        x = reshape(x,[],nSpecies);

        x = mirror(x);
        
        mu_s = GetChemicalPotential(x,t,mu);
        
        Dmu_s = Diff.Dy*mu_s;
        
        h_s   = Diff.Dy*x - Vext_grad;
        
        dxdt  = kBT*Diff.DDy*mu_s + h_s.*Dmu_s + 2*kBT*yS.^(-1).*Dmu_s; 
        
        dxdt(Ind.bound,:) = 0;

        dxdt = D0.*dxdt;
        
        dxdt = cut(dxdt);
        
        dxdt = dxdt(:);
 
    end
    
    function y = f(x)
        %solves for T*log*rho + Vext                
        mu_s         = x(1,:);
        mu_s         = repmat(mu_s,N,1);
        x            = x(2:end,:);
        x            = mirror(x);
        rho_full     = exp((x-Vext)/kBT);
        y            = GetChemicalPotential(x,0,mu_s);
        y            = cut(y);
        y            = [IntNP*rho_full - nParticlesS.';y];
    end
    
    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT); 
        mu_s = GetChemicalPotential(x,t,mu);
        flux  = -rho_s.*(Diff.Dy*mu_s);                                
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
        normalization = repmat( IntNP*rhoInit./nParticlesS' , size(rhoInit,1) ,1);
        x0=-kBT*log(normalization) - VAdd;
    end

    function w = NPweight(y)
        w = 4*pi*y.^2;
    end

    function mu = zeroFunction(rho_s,~,~)
        mu = zeros(size(rho_s));
    end

    function xOut = mirror(x)
        % flip for negative spacial part
        xOut = [flipdim(x,1); x];
    end

    function xOut = cut(x)
        % remove negative spacial part
        xOut = x(end/2+1:end,:);
    end


end