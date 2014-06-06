function [data,optsPhys,optsNum,optsPlot] = DDFT_Diffusion_1D_Planar(optsPhys,optsNum,optsPlot)

    if(nargin == 0)
        [data,optsPhys,optsNum,optsPlot] = Test_DDFT_DiffusionInfInterval();
        return;
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
    D0 = 1./(gamma.*m);
    
    plotTimes   = optsNum.plotTimes;
    
    Fex_Num = optsNum.FexNum;

    shapeClass = str2func(optsNum.PhysArea.shape);
    aLine = shapeClass(PhysArea);
    
    [Pts,Diff,Int,Ind,~] = aLine.ComputeAll(optsNum.PlotArea); 

    yS               = repmat(Pts.y,1,nSpecies);

    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************ 
        
    tic
    fprintf(1,'Computing Fex matrices ...');
    paramsFex.V2       = optsPhys.V2;
    paramsFex.kBT      = optsPhys.kBT;
    paramsFex.FexNum   = optsNum.FexNum;
    paramsFex.Pts      = aLine.Pts;     
    paramsFex.nSpecies = nSpecies;   
    
    if(strcmp(Fex_Num.Fex,'Percus'))
        IntMatrFex     = DataStorage(['FexData' filesep class(aLine)],@FexMatrices_Percus,paramsFex,aLine);
        getFex         = @Fex_Percus;
        getConv        = @zeroFunction;
        convStruct     = [];
    elseif(strcmp(Fex_Num.Fex,'Meanfield'))
        getFex         = @zeroFunction;
        IntMatrFex     = [];
        convStruct     = DataStorage(['FexData' filesep class(aLine)],@FexMatrices_Meanfield,paramsFex,aLine);
        getConv        = @Fex_Meanfield;
    elseif(strcmp(Fex_Num.Fex,'Perturbation'))
        IntMatrFex     = DataStorage(['FexData' filesep class(aLine)],@FexMatrices_Percus,paramsFex,aLine);
        getFex         = @Fex_Percus;
        convStruct     = DataStorage(['FexData' filesep class(aLine)],@FexMatrices_Meanfield,paramsFex,aLine);
        getConv        = @Fex_Meanfield;
    end
    fprintf(1,'done.\n');
    t_fex = toc;
    disp(['Fex computation time (sec): ', num2str(t_fex)]);
    
     if(isfield(optsNum,'HINum') && ~isempty(optsNum.HINum))
        doHI = true;
    else
        doHI = false;
    end
    
    if(doHI)
        tic
        fprintf(1,'Computing HI matrices ...');   
        paramsHI.optsPhys.HI       = optsPhys.HI;
        paramsHI.optsNum.HINum     = optsNum.HINum;
        paramsHI.optsNum.Pts       = aLine.Pts;
        paramsHI.optsPhys.nSpecies = nSpecies;
        HIStruct = DataStorage(['HIData' filesep class(aLine)],@HIMatrices,paramsHI,aLine);
        fprintf(1,'done.\n');
        t_HI = toc;
        display(['HI computation time (sec): ', num2str(t_HI)]); 
    else
        HIStruct = [];
    end
          
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    
    tic
    [Vext,Vext_grad]  = getVBackDVBack1D(yS,0,optsPhys.V1);
    
    x0 = getInitialGuess(0);   
    mu0 = zeros(1,nSpecies);
    x0mu0 = [mu0;x0];

    paramsIC.optsPhys.V1          = optsPhys.V1;
    paramsIC.optsPhys.V2          = optsPhys.V2;
    paramsIC.optsPhys.mS          = optsPhys.mS;
    paramsIC.optsPhys.kBT         = optsPhys.kBT;
    paramsIC.optsPhys.nParticlesS = optsPhys.nParticlesS;

    paramsIC.optsNum.FexNum   = optsNum.FexNum;
    paramsIC.optsNum.PhysArea = optsNum.PhysArea;
    
    fsolveOpts=optimset('Display','off');
    paramsIC.fsolveOpts = fsolveOpts;
    
    fprintf(1,'Computing initial condition ...');        
    eqStruct = DataStorage(['EquilibriumData' filesep class(aLine)],@ComputeEquilibrium,paramsIC,x0mu0);
    fprintf(1,'done.\n');
    
    x_ic = eqStruct.x_ic;
    flag = eqStruct.flag;
    
    if(flag<0)
        fprintf(1,'fsolve failed to converge\n');
        pause
    else
        fprintf(1,'Found initial equilibrium\n');
    end
    
    mu      = x_ic(1,:);
    mu      = repmat(mu,N,1);
    x_ic    = x_ic(2:end,:);
    
    t_eqSol = toc;
    disp(['Equilibrium computation time (sec): ', num2str(t_eqSol)]);
        
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    
    tic
    mM            = ones(N,1);
    mM(Ind.bound) = 0;
    mM = repmat(mM,nSpecies,1);
    opts = odeset('RelTol',10^-6,'AbsTol',10^-6,'Mass',diag(mM));       

    fprintf(1,'Computing dynamics ...'); 
    [~,X_t] =  ode15s(@dx_dt,plotTimes,x_ic,opts);
    
    fprintf(1,'done.\n');
    
    t_dynSol = toc;
    
    disp(['Dynamics computation time (sec): ', num2str(t_dynSol)]);

    %****************************************************************
    %****************  Post process                      ************
    %****************************************************************
    
    nPlots = length(plotTimes);
    
    X_t = X_t.';
    X_t = reshape(X_t,[],nSpecies,size(X_t,2));
    
    rho_t     = zeros(N,nSpecies,nPlots);
    flux_t    = zeros(N,nSpecies,nPlots);
    V_t       = zeros(N,nSpecies,nPlots);
    for i = 1:length(plotTimes)
        rho_t(:,:,i)  = exp((X_t(:,:,i)-Vext)/kBT);
        if(doHI)
            flux_t(:,:,i) = GetFlux_HI(X_t(:,:,i),plotTimes(i));
        else
            flux_t(:,:,i) = GetFlux(X_t(:,:,i),plotTimes(i));
        end

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
    function dxdt = dx_dt(t,x)
        
        x = reshape(x,[],nSpecies);
        
        mu_s  = GetChemicalPotential(x,t,mu);
        Dmu_s = Diff.Dy*mu_s;
        h_s   = Diff.Dy*x - Vext_grad;
        
        dxdt  = kBT*Diff.DDy*mu_s + (h_s.*Dmu_s); 
        
        if(doHI)
            rho_s = exp((x-Vext)/kBT);
            HI = ComputeHI(rho_s,Dmu_s,HIStruct);
            dxdt  = dxdt + kBT*Diff.Dy*HI + HI.*h_s;
        end
        
        flux_dir        = Diff.Dy*mu_s;
        
        if(doHI)
            flux_dir = flux_dir + HI;
        end
        
        dxdt(Ind.finite,:)     = Ind.normalFinite*flux_dir;
        dxdt(Ind.infinite,:)   = x(Ind.infinite,:) - x_ic(Ind.infinite,:);

        dxdt = D0.*dxdt;
        dxdt = dxdt(:);
 
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
    
    function flux = GetFlux_HI(x,t)
        rho_s = exp((x-Vext)/kBT); 
        mu_s = GetChemicalPotential(x,t,mu);
        Dmu_s = Diff.Dy*mu_s;
        HI = ComputeHI(rho_s,Dmu_s,HIStruct);
        flux  = -rho_s.*(Dmu_s + HI);                                
    end

    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT); 
        mu_s = GetChemicalPotential(x,t,mu);
        Dmu_s = Diff.Dy*mu_s;
        flux  = -rho_s.*Dmu_s;                                
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

    function eqStruct = ComputeEquilibrium(params,y0)      
        [y,status]   = fsolve(@f,y0,params.fsolveOpts);
        eqStruct.x_ic = y;
        eqStruct.flag = status;
    end

end