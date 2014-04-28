function [data,optsPhys,optsNum,optsPlot] = NS_SphericalInfInterval_old(optsPhys,optsNum,optsPlot)

    if(nargin == 0)
        [optsNum,optsPhys,optsPlot] = Default_NS_SphericalInfInterval();
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
    
    sigma = repmat(optsPhys.sigmaS',N,1);
    
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
    
    v_ic    = zeros(size(x_ic));
    xv_ic   = [x_ic ; v_ic];
   
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    mMx = ones(N/2,1);
    mMv = ones(N/2,1);
    mM  = [mMx;mMv]; 
    mM  = repmat(mM,nSpecies,1);
    opts = odeset('RelTol',10^-6,'AbsTol',10^-6,'Mass',diag(mM));
    [~,XV_t] =  ode15s(@dxv_dt,plotTimes,xv_ic,opts);        

    XV_t = XV_t.';
    XV_t = reshape(XV_t,[],nSpecies,size(XV_t,2));
    
    X_t  = XV_t(1:end/2,:,:);
    
    rho_t     = zeros(N,nSpecies,length(plotTimes));
    flux_t    = mirrorV(XV_t(end/2+1:end,:,:));
    V_t       = zeros(N,nSpecies,length(plotTimes));
    for i = 1:length(plotTimes)
        x_t = mirrorX(X_t(:,:,i));
        rho_t(:,:,i)  = exp((x_t-Vext)/kBT);
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
        
           
        NS = getNSterms(x,v);
        dvdt = dvdt - NS;  % NS terms calculated as if they're on the left of the equation

%         % Archer test
%         NSA = Diff.DDy*v + yS.^(-1).*(Diff.Dy*v) - 2*yS.^(-2).*v;
%         rho = exp((x - Vext)/kBT);
%         dvdt = dvdt + rho.^(2/3).*NSA;
        
        dxdt(Ind.bound,:) = 0;
        dvdt(Ind.bound,:) = 0;
       
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

    function NS = getNSterms(x,v)
        k0 = 1/6*(m*kBT).^(1/2)*pi.*sigma.^4;
        k1 = 1/30*(m*kBT).^(1/2)*pi.*sigma.^4;
        k2 = 1/5*(kBT*mInv).^(1/2)*pi;
        
        c1 = -k0 + 2/3*k1;
        c2 = 2/3*k2;
        %d1 = -2*k1;
        %d2 = -2*k2;

        d1 = 0;
        d2 = 0;
        
        Dv  = Diff.Dy*v;
        DDv = Diff.DDy*v;
        
                        
        DlnRho = (Diff.Dy*x - Vext_grad)/kBT;
       
        rho = exp((x-Vext)/kBT);
        rho23 = rho.^(2/3);
        
        NS =  (c1.*rho23 + c2).*(2*yS.^(-2) .*v + 4 * yS.^(-1) .* Dv + DDv) ...
            + (d1.*rho23 + d2).*(DDv + 2 * yS.^(-1) .* Dv) ...
             + (5/3*c1.*rho23 + c2).*DlnRho.*( 2*yS.^(-1).*v  + Dv) ...
             + (5/3*d1.*rho23 + d2).*DlnRho.*Dv;
        
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
        % flip for negative spatial part and change sign
        vOut = [-flipdim(v,1); v];
    end

    function xOut = cut(x)
        % remove negative spatial part
        xOut = x(end/2+1:end,:);
    end


end