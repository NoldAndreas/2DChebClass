function data = DDFT_DiffusionPlanar_NSpecies_Andreas(optsPhys,optsNum,optsPlot)
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
    if(nargin==0)
        %[optsNum,optsPhys] = Test_DDFT_DiffusionPlanar_NSpecies();
        [optsNum,optsPhys,optsPlot] = TestFMT_DDFT_DiffusionPolar_NSpecies(false);
    end

    disp(['** ',optsNum.DDFTCode,' **']);
    
    %************************************************
    %***************  Initialization ****************
    %************************************************
        
    [kBT,nParticlesS,HS_f]           = LoadPhysData_2Species(optsPhys);
    [N1,N2,PhysArea,plotTimes]       = LoadNumDataLoc(optsNum);    
    optsPhys.FexNum  = optsNum.FexNum;
              
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    
    tic                  
    infSpace                  = InfSpace(PhysArea);
    [Pts,Diff,Int,Ind,Interp] = infSpace.ComputeAll(optsNum.PlotArea);                                        
	convStruct     = DataStorage('FexData',@FexMatrices_Meanfield,optsPhys,infSpace);   
    %convStruct                = FexMatrices_Meanfield(optsPhys,infSpace);          
    
    infSpace.doPlots(convStruct(1,1).Conv*ones(size(Pts.y1_kv)));
    %doPlots_SC(Interp,Pts,convStruct(1,1).Conv*ones(size(Pts.y1_kv)));
     
    nSpecies=length(nParticlesS);
    
    y1S               = repmat(Pts.y1_kv,1,nSpecies); 
    y2S               = repmat(Pts.y2_kv,1,nSpecies);
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
     fsolveOpts=optimset('Display','off');
     [x_ic,flag]    = fsolve(@f,y0,fsolveOpts); 
     
     if(flag<0)
         fprintf(1,'fsolve failed to converge\n');
         pause
     else
         fprintf(1,'Found initial equilibrium\n');
     end
      
%     x_ic    = fsolve(@f,zeros(1+N1*N2,nSpecies)); 
    
    mu      = x_ic(1,:);
    x_ic    = x_ic(2:end,:);        
    
    rho_ic  = exp((x_ic-Vext)/kBT);       
    
    if(nargin<3)
        optsPlot.lineColourDDFT={'r','b','g','k','m'};
        optsPlot.doDDFTPlots=true;
    end
    
    if(optsPlot.doDDFTPlots)
        figure
        infSpace.doPlots(rho_ic,'',optsPlot.lineColourDDFT);
        %doPlots_IP_Contour_NSpecies(Interp,rho_ic,optsPlot)
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
    
    [outTimes,X_t] =  ode15s(@dx_dt,plotTimes,x_ic,opts);        
    t_dynSol = toc;

    %************************************************
    %****************  Postprocess  ****************
    %************************************************        

    X_t       = X_t';
        
    rho_t     = exp((X_t-Vext(:)*ones(1,length(outTimes)))/kBT);
    
    X_t       = reshape(X_t,N1*N2,nSpecies,length(outTimes));
    rho_t     = reshape(rho_t,N1*N2,nSpecies,length(outTimes));
    flux_t    = zeros(2*N1*N2,nSpecies,length(plotTimes));
    
%     for i = 1:length(plotTimes)
%         flux_t(:,:,i) = GetFlux(X_t(:,:,i),plotTimes(i));        
%     end
    
    data       = v2struct(convStruct,X_t,rho_t,...
                        mu,flux_t,...
                        t_preprocess,t_eqSol,t_dynSol);                        
    data.shape = infSpace;
    
    if(~isfield(optsNum,'savefileDDFT'))
    %    SaveToFile(optsNum.DDFTCode,data,optsPhys,optsNum,getResultsPath());
    end                    
                    
	display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
    display(['Equilibrium Sol., Computation time (sec): ', num2str(t_eqSol)]);
    display(['Dynamics, Computation time (sec): ', num2str(t_dynSol)]);
     
    if(optsPlot.doDDFTPlots) 
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
        mu_s(Pts.y1_kv==inf | Pts.y2_kv==inf,:) = 0;
        h_s      = Diff.grad*x - Vext_grad;
        
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        
        dxdt(Ind.bound,:)  = x(Ind.bound,:) - x_ic(Ind.bound,:);
        
        %Boundary Conditions: no flux at the walls        
        %flux_dir          = Diff.grad*mu_s;
        %dxdt(Ind.bound,:) = Ind.normal*flux_dir;           
             
        dxdt = dxdt(:);
    end

    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp((x-Vext)/kBT);
                 
        mu_s = Fex_Meanfield(rho_s,convStruct,kBT);%getFex(rho_s,IntMatrFex,kBT);
                       
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

    %***************************************************************
    %Auxiliary functions:
    %***************************************************************
    function [kBT,nParticles,HS_f] = LoadPhysData_2Species(optsPhys)        
        kBT        = optsPhys.kBT;   
        nParticles = optsPhys.nParticlesS;     
        if(isfield(optsPhys,'HSBulk'))
            HS_f       = str2func(optsPhys.HSBulk);  
        elseif(isfield(optsNum,'HSBulk'))
            HS_f       = str2func(optsNum.HSBulk); 
        else
            HS_f       = @ZeroMap;
        end     
    end

    function [N1,N2,PhysArea,plotTimes] = LoadNumDataLoc(optsNum)
        N1 = optsNum.PhysArea.N1;
        N2 = optsNum.PhysArea.N2;
        PhysArea = optsNum.PhysArea;
        PhysArea.N = [N1,N2];
        PhysArea.L1 = optsNum.PhysArea.L1;
        PhysArea.L2 = optsNum.PhysArea.L2;
        plotTimes  = optsNum.plotTimes;
    end

    function [muSC,fnCS] = ZeroMap(rho,kBT)
        muSC = 0;
        fnCS = 0;
    end

    

end