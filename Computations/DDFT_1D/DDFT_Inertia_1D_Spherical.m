function [data,optsPhys,optsNum,optsPlot] = DDFT_Inertia_1D_Spherical(optsPhys,optsNum,optsPlot)

    if(nargin == 0)
        [data,optsNum,optsPhys,optsPlot] = Test_DDFT_InertiaSphericalInfInterval();
        return;
    end
    
    AddPaths();
    
    %close all;  
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
    
    mInv  = m.^(-1);
    
    plotTimes   = optsNum.plotTimes;
    
    Fex_Num = optsNum.FexNum;
    
    shapeClass = str2func(optsNum.PhysArea.shape);
    aLine = shapeClass(PhysArea);
    
    [Pts,Diff,~,Ind,~] = aLine.ComputeAll(optsNum.PlotArea); 

    yS               = repmat(Pts.y,1,nSpecies);
            
    %************************************************
    %****************  Preprocess  ****************
    %************************************************ 
    
    % SORT THE HARD CODING HERE
    geom.yMin = 0;
    geom.yMax = 20;
    IntNP = aLine.IntegrateRegion(@NPweight,geom);
       
    yS               = repmat(Pts.y,1,nSpecies);
    
    tic
    fprintf(1,'Computing Fex matrices ...');
    paramsFex.V2       = optsPhys.V2;
    paramsFex.kBT      = optsPhys.kBT;
    paramsFex.FexNum   = optsNum.FexNum;
    paramsFex.Pts      = aLine.Pts;     
    paramsFex.nSpecies = nSpecies;   
    
    if(strcmp(Fex_Num.Fex,'FMT'))
        IntMatrFex     = DataStorage(['FexData' filesep class(aLine)],@FexMatrices_FMTRosenfeldSpherical,paramsFex,aLine);
        getFex         = @Fex_FMTRosenfeldSpherical;
        getConv        = @zeroFunction;
        convStruct     = [];
    elseif(strcmp(Fex_Num.Fex,'Meanfield'))
        getFex         = @zeroFunction;
        IntMatrFex     = [];
        convStruct     = DataStorage(['FexData' filesep class(aLine)],@FexMatrices_Meanfield_Spherical,paramsFex,aLine);
        getConv        = @Fex_Meanfield;
    elseif(strcmp(Fex_Num.Fex,'Perturbation'))
        IntMatrFex     = DataStorage(['FexData' filesep class(aLine)],@FexMatrices_FMTRosenfeldSpherical,paramsFex,aLine);
        getFex         = @Fex_FMTRosenfeldSpherical;
        convStruct     = DataStorage(['FexData' filesep class(aLine)],@FexMatrices_Meanfield_Spherical,paramsFex,aLine);
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
        HIStruct = DataStorage(['HIData' filesep class(aLine)],@HIMatricesSpherical,paramsHI,aLine);
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
    
    % add in VAdd here
    x0 = getInitialGuess(0);   
    x0 = cut(x0);
    
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
    
    %x_ic    = fsolve(@f,x0mu0);
    mu      = x_ic(1,:);
    mu      = repmat(mu,N,1);
    x_ic    = x_ic(2:end,:);
    x_icFull = mirrorX(x_ic);
    
    v_ic    = zeros(size(x_ic));
    xv_ic   = [x_ic ; v_ic];
    
    t_eqSol = toc;
    disp(['Equilibrium computation time (sec): ', num2str(t_eqSol)]);
   
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    mMx = ones(N/2,1);
    mMv = ones(N/2,1);
    mMx(Ind.infinite(N/2+1:end)) = 0;
    mMv(Ind.bound(N/2+1:end)) = 0;
    mM  = [mMx;mMv]; 
    mM  = repmat(mM,nSpecies,1);
    opts = odeset('RelTol',10^-6,'AbsTol',10^-6,'Mass',diag(mM));
    %[~,XV_t] =  ode15s(@dxv_dt,plotTimes,xv_ic,opts);        

    fprintf(1,'Computing dynamics ...'); 
    [~,XV_t] =  ode15s(@dxv_dt,plotTimes,xv_ic,opts);
    fprintf(1,'done.\n');
    
    t_dynSol = toc;
    
    disp(['Dynamics computation time (sec): ', num2str(t_dynSol)]);
    
    %****************************************************************
    %****************  Post process                      ************
    %****************************************************************
    
    nPlots = length(plotTimes);
    
    XV_t = XV_t.';
    XV_t = reshape(XV_t,[],nSpecies,size(XV_t,2));
    
    X_t  = XV_t(1:end/2,:,:);
    
    rho_t     = zeros(N,nSpecies,nPlots);
    flux_t    = mirrorV(XV_t(end/2+1:end,:,:));
    V_t       = zeros(N,nSpecies,nPlots);
    for i = 1:length(plotTimes)
        x_t = mirrorX(X_t(:,:,i));
        rho_t(:,:,i)  = exp((x_t-Vext)/kBT);
        V_t(:,:,i)    = Vext + getVAddDVAdd1D(yS,plotTimes(i),optsPhys.V1);
    end
    
    data       = v2struct(IntMatrFex,convStruct,HIStruct,X_t,rho_t,mu,flux_t,V_t);
    data.shape = aLine;
    data.shape.Int = IntNP;
    
    if(~isfield(optsNum,'doPlots') ...
            || (isfield(optsNum,'doPlots') && optsNum.doPlots) )
        PlotDDFT1D(v2struct(optsPhys,optsNum,data));       
    end
    
   %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function dxvdt = dxv_dt(t,xv)
        
        xv = reshape(xv,[],nSpecies);

        x = xv(1:end/2,:);
        v = xv(end/2+1:end,:);
        
        x = mirrorX(x);
        v = mirrorV(v);
        
        mu_s = GetChemicalPotential(x,t,mu);
        
        h_s   = Diff.Dy*x - Vext_grad;
        
        dxdt  = - kBT*Diff.Dy*v - h_s.*v - 2*kBT*yS.^(-1).*v;  
        
        dvdt  = - gamma.*v - v.*(Diff.Dy*v) - mInv.*(Diff.Dy*mu_s);
        
        if(doHI)
            rho_s = exp((x-Vext)/kBT);
            HI = ComputeHI(rho_s,v,HIStruct);
            dvdt = dvdt - gamma.*HI;
        end

        dxdt(Ind.infinite,:) = x(Ind.infinite,:) - x_icFull(Ind.infinite,:);
        dvdt(Ind.bound,:) = v(Ind.bound,:);

        dxdt = cut(dxdt);
        dvdt = cut(dvdt);
        
        dxvdt = [dxdt;dvdt];
        
        dxvdt = dxvdt(:);
 
    end

    function y = f(x)
        %solves for T*log*rho + Vext                
        mu_s         = x(1,:);
        mu_s         = repmat(mu_s,N,1);
        x            = x(2:end,:);
        x            = mirrorX(x);
        rho_full     = exp((x-Vext)/kBT);
        y            = GetChemicalPotential(x,0,mu_s);
        y            = cut(y);
        y            = [IntNP*rho_full - nParticlesS.';y];
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

    function xOut = mirrorX(x)
        % flip for negative spatial part
        xOut = [flipdim(x,1); x];
    end

    function vOut = mirrorV(v)
        % flip for negative velocity part and change sign
        vOut = [-flipdim(v,1); v];
    end

    function xOut = cut(x)
        % remove negative spatial part
        xOut = x(end/2+1:end,:);
    end

    function eqStruct = ComputeEquilibrium(params,y0)      
        [y,status]   = fsolve(@f,y0,params.fsolveOpts);
        eqStruct.x_ic = y;
        eqStruct.flag = status;
    end


end