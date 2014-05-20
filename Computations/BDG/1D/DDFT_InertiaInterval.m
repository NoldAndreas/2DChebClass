function [data,optsPhys,optsNum,optsPlot] = DDFT_InertiaInterval(optsPhys,optsNum,optsPlot)

    if(nargin == 0)
        [data,optsNum,optsPhys,optsPlot] = Test_DDFT_InertiaInterval();
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

    mInv  = m.^(-1);
    
    plotTimes   = optsNum.plotTimes;
    
    Fex_Num = optsNum.FexNum;
            
    %************************************************
    %****************  Preprocess  ****************
    %************************************************ 
    
    aLine = SpectralLine(PhysArea);

    [Pts,Diff,Int,Ind,~] = aLine.ComputeAll(optsNum.PlotArea); 

    yS               = repmat(Pts.y,1,nSpecies);
    
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
    [x_ic,flag]     = DataStorage(['EquilibriumData' filesep class(aLine)],@ComputeEquilibrium,paramsIC,x0mu0);
    fprintf(1,'done.\n');
    
    if(flag<0)
        fprintf(1,'fsolve failed to converge\n');
        pause
    else
        fprintf(1,'Found initial equilibrium\n');
    end

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
    mMv(Ind.bound) = 0;
    mM             = [mMx;mMv];
    mM = repmat(mM,nSpecies,1);
    opts = odeset('RelTol',10^-4,'AbsTol',10^-4,'Mass',diag(mM));
    
    fprintf(1,'Computing dynamics ...'); 
    [~,XV_t] =  ode15s(@dxv_dt,plotTimes,xv_ic,opts);
    fprintf(1,'done.\n');
    t_dynSol = toc;
    disp(['Dynamics computation time (sec): ', num2str(t_dynSol)]);

    %****************************************************************
    %****************  Post process                      ************
    %****************************************************************


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
    
    data       = v2struct(IntMatrFex,convStruct,X_t,rho_t,mu,flux_t,V_t);
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
        
        %Boundary Conditions: no flux at the walls
        dvdt(Ind.bound,:) = v(Ind.bound,:);

        dxvdt = [dxdt;dvdt];
        
        dxvdt = dxvdt(:);         
 
    end
    
    function y = f(x)
        %solves for T*log*rho + Vext                
        mu_s         = x(1,:);
        mu_s         = repmat(mu_s,N,1);
        x            = x(2:end,:);        
        rho_full     = exp((x-Vext)/kBT);
        y            = GetChemicalPotential(x,0,mu_s); %y = x + Conv*rho_full - mu_s + getVAdd(Pts.y1_kv,Pts.y2_kv,0,optsPhys);                
        y            = [Int*rho_full - nParticlesS.';y];
    end

    function mu_s = GetChemicalPotential(x,t,mu)
        rho_s = exp((x-Vext)/kBT); 
        mu_s  = getFex(rho_s,IntMatrFex,kBT);
        VAdd  = getVAddDVAdd1D(yS,t,optsPhys.V1);
        %ConvRho_s = Fex_Meanfield(rho_s,convStruct,kBT);
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

    function [y,flag] = ComputeEquilibrium(params,y0)      
        [y,flag]   = fsolve(@f,y0,params.fsolveOpts); 
    end

end