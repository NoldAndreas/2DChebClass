function data = DDFT_DiffusionPlanar_NSpecies(optsPhys,optsNum,optsPlot)
%************************************************************************* 
% data = DDFT_DiffusionPlanar_NSpecies(optsPhys,optsNum,optsPlot)
%
% define
%   mu_s_i := kBT*log(rho_i) + sum( int(rho_j(r')*Phi2D(r-r'),dr') , 
%               j = 1..N) +mu_HS(rho_i) + V_ext_i - mu_i
% Equilibrium, (solve for each species i)
%  (EQ i,1) mu_s_i     = 0
%  (EQ i,2) int(rho_i) = NoParticles_i
% Dynamics:
%   (DYN i) rho_i/dt = div(rho_i*grad(mu_s_i))
%*************************************************************************   
    if(nargin == 0)
        [optsNum,optsPhys,optsPlot] = TestFMT_DDFT_DiffusionPolar_NSpecies(false);
    end

    disp(['** ',optsNum.DDFTCode,' **']);

    %************************************************
    %***************  Initialization ****************
    %************************************************
    PhysArea    = optsNum.PhysArea;   
    N1 = PhysArea.N(1); N2 = PhysArea.N(2);              
    kBT          = optsPhys.kBT; 
    nParticlesS  = optsPhys.nParticlesS;            
    plotTimes    = optsNum.plotTimes;
    
    if(isfield(optsPhys,'HSBulk'))
        HS_f       = str2func(optsPhys.HSBulk);  
    elseif(isfield(optsNum,'HSBulk'))
        HS_f       = str2func(optsNum.HSBulk); 
    else
        HS_f       = @ZeroMap;
    end                     
    getFex = str2func(['Fex_',optsNum.FexNum.Fex]);
    
    if(nargin<3)
        optsPlot.lineColourDDFT={'r','b','g','k','m'};
        optsPlot.doDDFTPlots=true;
    end
    
    L1   = PhysArea.L1; L2   = PhysArea.L2; N = [N1,N2];
    IDC = InfSpace(v2struct(L1,L2,N));
    
    [Pts,Diff,Int,Ind,Interp] = IDC.ComputeAll(optsNum.PlotArea);
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************           
    tic
    fprintf(1,'Computing Fex matrices ...');   
    params.V2          = optsPhys.V2;
    params.kBT         = optsPhys.kBT;
    params.FexNum      = optsNum.FexNum;
    params.Pts         = IDC.Pts;     
    params.Polar       = 'cart';
    params.nParticlesS = optsPhys.nParticlesS;
       
    func           = str2func(['FexMatrices_',optsNum.FexNum.Fex]);    
	IntMatrFex     = DataStorage('FexData',func,params,IDC);   
        
    fprintf(1,'done.\n');
    tfex = toc;
    disp(tfex);
    
    if(isfield(optsPhys,'doHI'))
        doHI = optsPhys.doHI;
    else
        doHI = false;
    end
    if(doHI)        
        tic
        fprintf(1,'Computing HI matrices ...');   
        paramsHI.kBT  = optsPhys.kBT;
        paramsHI.HI   = optsPhys.HI;
        paramsHI.HIShapeParams  = optsNum.HIShapeParams;
        paramsHI.Pts     = IDC.Pts;     
        paramsHI.Polar   = 'cart';
        func           = str2func('HIMatrices_RP');    
        IntMatrHI     = DataStorage('HIData',func,paramsHI,IDC);   

        %IntMatrHI = HIMatrices_RP(optsPhys,optsNum,IDC);

        fprintf(1,'done.\n');
        tHI = toc;
        disp(tHI);        
    end
    
    nSpecies=length(nParticlesS);
    
    y1S     = repmat(Pts.y1_kv,1,nSpecies); 
    y2S     = repmat(Pts.y2_kv,1,nSpecies);
    [Vext,Vext_grad]  = getVBackDVBack(y1S,y2S,optsPhys.V1);       

    I  = eye(N1*N2);
    eyes=repmat(I,1,2);
    
    t_preprocess = toc;
  
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    
    tic
    
    VAdd0=getVAdd(y1S,y2S,0,optsPhys.V1);
    y0 = getInitialGuess(VAdd0);
    
    %y0=zeros(1+N1*N2,nSpecies);
    
%     x_ic    = fsolve(@f,zeros(1+N1*N2,nSpecies));
    fprintf(1,'Computing initial condition ...');
    fsolveOpts=optimset('Display','off');
    [x_ic,flag]    = fsolve(@f,y0,fsolveOpts); 
    fprintf(1,'done.\n');
    
    if(flag<0)
        fprintf(1,'fsolve failed to converge\n');
        pause
    else
        fprintf(1,'Found initial equilibrium\n');
    end
    
    mu      = x_ic(1,:);
    x_ic    = x_ic(2:end,:);
    
    if(optsPlot.doDDFTPlots)
        figure
        rho_ic  = exp((x_ic-Vext)/kBT);
        IDC.doPlots(rho_ic,'',optsPlot.lineColourDDFT);        
    end

    t_eqSol = toc;
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
    
    tic
    mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    mM=repmat(mM,nSpecies,1);
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    
    if(doHI)
        [outTimes,X_t] =  ode15s(@dx_dt_HI,plotTimes,x_ic,opts);        
    else
        [outTimes,X_t] =  ode15s(@dx_dt,plotTimes,x_ic,opts);
    end
    t_dynSol = toc;
    
    %************************************************
    %****************  Postprocess  ****************
    %************************************************        

    X_t       = X_t';
    rho_t     = exp((X_t-Vext(:)*ones(1,length(outTimes)))/kBT);
    
    X_t       = reshape(X_t,N1*N2,nSpecies,length(outTimes));
    rho_t     = reshape(rho_t,N1*N2,nSpecies,length(outTimes));
    flux_t    = zeros(2*N1*N2,nSpecies,length(plotTimes));
    
    for i = 1:length(plotTimes)
        if(doHI)
            flux_t(:,:,i) = GetFlux_HI(X_t(:,:,i),plotTimes(i));        
        else
            flux_t(:,:,i) = GetFlux(X_t(:,:,i),plotTimes(i));        
        end
    end
    
    data     = v2struct(IntMatrFex,X_t,rho_t,mu,flux_t,...
                        t_preprocess,t_eqSol,t_dynSol);                        
    data.shape = IDC;
    
    %if(~isfield(optsNum,'savefileDDFT'))
    %    SaveToFile(optsNum.DDFTCode,data,optsPhys,optsNum,getResultsPath());
    %end                    
                    
	display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
    display(['Equilibrium Sol., Computation time (sec): ', num2str(t_eqSol)]);
    display(['Dynamics, Computation time (sec): ', num2str(t_dynSol)]);
    
    
    if(~exist('optsPlot','var') || optsPlot.doDDFTPlots)
        PlotDDFT(v2struct(optsPhys,optsNum,optsPlot,data));  
    end
             
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function y = f(x)
        %solves for T*log*rho + Vext                
        mu_s         = x(1,:);
        x            = x(2:end,:);       
        rho_full     = exp((x-Vext)/kBT);
        y            = GetExcessChemPotential(x,0,mu_s); 
        y            = [Int*rho_full - nParticlesS';y];
        y            = y(:);
    end

    function y0 = getInitialGuess(VAdd)
        
      if(~isfield(optsPhys,'HSBulk') && ~isfield(optsNum,'HSBulk'))
        % initial guess for mu doesn't really matter
        muInit=zeros(1,nSpecies);

        % however, require relatively good guess for rho

        % rho without interaction ie exp(-(VBack+VAdd )/kBT)
        rhoInit = exp(-(Vext+VAdd)/kBT);
        
        % inverse of normalization coefficient
        normalization = repmat( Int*rhoInit./nParticlesS' , size(rhoInit,1) ,1);

        %rho     = exp((x-Vext)/kBT) = N exp((-Vadd - Vext)/kBT
        %x/kBT - Vext/kBT = log(N) - Vadd/kBT - Vext/kBT
        %x = kBT log(N) - Vadd
        
        y0=-kBT*log(normalization) -VAdd;
        
        y0 = [muInit; y0];
      else
          y0 = zeros(1+N1*N2,nSpecies);
      end
    
    end

    function dxdt = dx_dt(t,x)
        x       = reshape(x,N1*N2,nSpecies);
        
        mu_s     = GetExcessChemPotential(x,t,mu);
        mu_s(Ind.bound,:) = 0;
        h_s      = Diff.grad*x - Vext_grad;
        h_s(Ind.bound,:) = 0; %here, we have assumed that grad(mu) converges fast enough
                
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        
        %Boundary Conditions at infinity
 %       if(PhysArea.y1Max == Inf)        
            dxdt(Ind.bound,:)  = x(Ind.bound,:) - x_ic(Ind.bound,:);   
        %else
%            flux_dir        = Diff.grad*mu_s;
%            dxdt(Ind.outR,:)  = Ind.normalOutR*flux_dir;             
%        end
        
        dxdt = dxdt(:);
    end

    function dxdt = dx_dt_HI(t,x)
        x       = reshape(x,N1*N2,nSpecies);
        
        mu_s     = GetExcessChemPotential(x,t,mu);
        mu_s(Ind.bound,:) = 0;
        
        rho_s = exp((x-Vext)/kBT);
        gradMu_s = Diff.grad*mu_s;
        HI_s = OverdampedHI(rho_s,IntMatrHI,gradMu_s);
        
        h_s      = Diff.grad*x - Vext_grad;
        h_s(Ind.bound,:) = 0; %here, we have assumed that grad(mu) converges fast enough
                
        dxdt     = kBT*(Diff.Lap*mu_s + Diff.div*HI_s) + eyes*( h_s.*(gradMu_s + HI_s) );  
          
        dxdt(Ind.bound,:)  = x(Ind.bound,:) - x_ic(Ind.bound,:);   

        
        dxdt = dxdt(:);
    end

    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp((x-Vext)/kBT);
        
        mu_s = getFex(rho_s,IntMatrFex,kBT);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end

        mu_s = mu_s + x + HS_f(rho_s,kBT) + getVAdd(y1S,y2S,t,optsPhys.V1);
                   
    end

    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = GetExcessChemPotential(x,t,mu); 
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end

    function flux = GetFlux_HI(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = GetExcessChemPotential(x,t,mu); 
        gradMu_s = Diff.grad*mu_s;
        HI_s = OverdampedHI(rho_s,IntMatrHI,gradMu_s);
        flux  = -[rho_s;rho_s].*(gradMu_s + HI_s);                                
    end


    %***************************************************************
    %Auxiliary functions:
    %***************************************************************

    function [muSC,fnCS] = ZeroMap(h1s,h2s)
        muSC = 0;
        fnCS = 0;
    end

    

end