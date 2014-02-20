function data = DDFT_DiffusionInfSpace_NSpecies(optsPhys,optsNum,optsPlot,name)
%************************************************************************* 
% define
%   mu_s_i := kBT*log(rho_i) + sum( int(rho_j(r')*Phi2D(r-r'),dr') , 
%               j = 1..N) +mu_HS(rho_i) + V_ext_i - mu_i
% Equilibrium, (solve for each species i)
%  (EQ i,1) mu_s_i     = 0
%  (EQ i,2) int(rho_i) = NoParticles_i
% Dynamics: 
%*************************************************************************   
    if(nargin == 0)
        Phys_Area = struct('N',[15;15],'L1',2,'L2',2.,'y2wall',0.,...
                           'N2bound',16,'h',1,'L2_AD',1.);

        Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                           'y2Min',-5,'y2Max',5,'N2',100);

        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                           'Ncircle',1,'N1disc',24,'N2disc',24);

        optsNum = struct('PhysArea',Phys_Area,...
                         'PlotArea',Plot_Area,...
                         'FexNum',Fex_Num,...
                         'DDFTCode','DDFT_DiffusionInfSpace_NSpecies',...
                         'Tmax',3,'TN',50,...
                         'name','default',...
                         'Accuracy_Averaging',1e-6);

        %V1       = struct('V1DV1','Vext_Cart_1','y10',0,'y20',0,'V0',2,'grav',1);              
        V1       = struct('V1DV1','Vext_Cart_6','V0',0,'L1',1,'L2',1,'Axis','y1y2','grav_start',0,'grav_end',-1,'T',1,'L1add',1,'L2add',1);
        %V1       = struct('V1DV1','zeroPotential');
        V2       = struct('V2DV2','hardSphere','sigmaS',1);                                 

        optsPhys = struct('V1',V1,'V2',V2,...                                            
                          'kBT',1,'eta',0.4,...%4783,...%0.4257,...
                          'nParticlesS',10);

        lineColourDDFT={{'r','b','g'}};            
        optsPlot = struct('lineColourDDFT',lineColourDDFT);
        optsPlot.doDDFTPlots=true;
    end

    close all;  
    disp(['** ',optsNum.DDFTCode,' **']);

    %************************************************
    %***************  Initialization ****************
    %************************************************        
    [kBT,nParticlesS,HS_f,eta]  = LoadPhysData_2Species(optsPhys);
    PhysArea                    = optsNum.PhysArea;   
    nSpecies                    = length(nParticlesS);
    plotTimes = 0:(optsNum.Tmax/optsNum.TN):optsNum.Tmax;
    R         = diag(optsPhys.V2.sigmaS)/2;
    fBulk     = str2func(['FexBulk_',optsNum.FexNum.Fex]);
    
    N1        = PhysArea.N(1);   N2 = PhysArea.N(2);
        
    getFex = str2func(['Fex_',optsNum.FexNum.Fex]);
    
    if(nargin<3)
        optsPlot.lineColourDDFT={'r','b','g','k','m'};
        optsPlot.doDDFTPlots=true;
    end
    
    
    IS                        = InfSpace_FMT(PhysArea);
    [Pts,Diff,Int,Ind,Interp] = IS.ComputeAll(optsNum.PlotArea);    
    
    rhoBulk   = eta*6/pi;
    mu        = kBT*log(rhoBulk) + fBulk(rhoBulk,kBT);
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************
    tic
    fprintf(1,'Computing Fex matrices ...');   
    
    params         = optsPhys;  
    params         = rmfield(params,'V1');
    params.FexNum  = optsNum.FexNum;
    params.PhysArea = optsNum.PhysArea;
    params.Pts     = IS.Pts;     
    params.Polar   = 'cart';      
    func           = str2func(['FexMatrices_',optsNum.FexNum.Fex]);       
    IntMatrFex     = DataStorage([class(IS) filesep func2str(func)],func,params,IS); %true    
        
    fprintf(1,'done.\n');
    tfex = toc;
    disp(tfex);
    
    y1S     = repmat(Pts.y1_kv,1,nSpecies); 
    y2S     = repmat(Pts.y2_kv,1,nSpecies);
    [Vext,Vext_grad]  = getVBackDVBack(y1S,y2S,optsPhys.V1);       

    I       = eye(N1*N2);
    eyes    = repmat(I,1,2);
    
    t_preprocess = toc;
    
    %*****************************************************
    %**************** Plot external potential   **********
    %*****************************************************
%      for k = 1:length(plotTimes)
%          V =  getVAdd(y1S,y2S,plotTimes(k),optsPhys.V1);%Vext +
%          HS.doPlots(V,true,false);
%          title(['t = ',num2str(plotTimes(k))]);
%          xlim([-5 5]);
%          ylim([0.5 5]);
%          pause(0.02);        
%      end
%   
    %****************************************************************
    %**************** Solve for equilibrium 2D condition   **********
    %****************************************************************
    disp('Compute 3D profile..');
        
    markComp       = true(size(Pts.y1_kv));    
    x_ic           = zeros(size(Pts.y1_kv));    
    icParams       = v2struct(optsPhys,optsNum);
    %icParams.optsPhys.V1 = rmfield(icParams.optsPhys.V1,'grav_end');
    %icParams.optsNum     = rmfield(icParams.optsNum,'Tmax');
    %icParams.optsNum     = rmfield(icParams.optsNum,'TN');
    
    icParams.fsolveOpts = optimset();%'MaxFunEvals',2000000,'MaxIter',200000);    
    x_ic           = DataStorage('EquilibriumData',@ComputeEquilibrium,icParams,x_ic);         	
    rho_ic         = exp((x_ic-Vext)/kBT);
        
    figure;    
	IS.doPlots(rho_ic,'SC');
	%PlotRosenfeldFMT_AverageDensities(IS,IntMatrFex(1),rho_ic);        
    
    %****************************************************************
    %**************** Solve for 3D dynamics *************************
    %****************************************************************
    mM            = ones(N1*N2,1);
    %mM(Ind.bound) = 0;
    opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
    [outTimes,X_t] =  ode15s(@dx_dt,[0 max(plotTimes)],x_ic,opts);                
        
    plotTimes = outTimes;
    X_t       = X_t';
    rho_t     = exp((X_t-Vext*ones(1,length(outTimes)))/kBT);
    
    flux_t    = zeros(2*N1*N2,length(plotTimes));
    for i = 1:length(plotTimes)
        flux_t(:,i) = GetFlux(X_t(:,i),plotTimes(i));
    end
    
    data              = v2struct(X_t,rho_t,mu,flux_t);%Pts,Diff,Int,Ind,Interp
    data.shape        = IS;
    optsNum.plotTimes = plotTimes;  
    PlotDDFT(v2struct(optsPhys,optsNum,data));                 

    %********************************************
    %**************** Postprocessing   **********
    %********************************************
    t_eqSol = toc;

	display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
    display(['Equilibrium Sol., Computation time (sec): ', num2str(t_eqSol)]);

    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************                  
    function y = ComputeEquilibrium(params,y0)
        fprintf(1,'Computing initial condition ...');        
        y          = fsolve(@f,y0,params.fsolveOpts);     
    end
    function y = f(x)
        %solves for T*log*rho + Vext                        
        y            = GetExcessChemPotential(x,0,mu);         
        y            = y(:);
    end

    function mu_s = GetExcessChemPotential(x,t,mu_offset)
        rho_s = exp((x-Vext(markComp))/kBT);                
        mu_s  = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - mu_offset(iSpecies);
        end
        Vadd = getVAdd(y1S,y2S,t,optsPhys.V1);
        mu_s = mu_s + x  + Vadd(markComp);  %HS_f(rho_s,kBT)
    end

    function dxdt = dx_dt(t,x)
       % x       = reshape(x,N1*N2,nSpecies);
        
        mu_s     = GetExcessChemPotential(x,t,mu);
        mu_s((Pts.y1_kv==inf) | (Pts.y1_kv==-inf) | (Pts.y2_kv==inf) | (Pts.y2_kv==-inf),:) = 0;
        h_s      = Diff.grad*x - Vext_grad;
        %h_s([Pts.y2_kv==inf;Pts.y2_kv==inf],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([false(N1*N2,1);Pts.y2_kv==inf],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([false(N1*N2,1);Pts.y2_kv==-inf],:) = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==inf;false(N1*N2,1)],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==-inf;false(N1*N2,1)],:) = 0; %here, we have assumed that grad(mu) converges fast enough
                
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
                
        %Boundary Conditions at infinity
        dxdt(Ind.bound,:)  = x(Ind.bound,:) - x_ic(Ind.bound,:);   
        
        flux_dir           = Diff.grad*mu_s;                
        dxdt(Ind.bottom,:) = Ind.normalBottom*flux_dir;             
               
        dxdt = dxdt(:);      
    end

    %***************************************************************
    %Auxiliary functions:
    %***************************************************************
    function [kBT,nParticles,HS_f,eta] = LoadPhysData_2Species(optsPhys)        
        kBT        = optsPhys.kBT;   
        nParticles = optsPhys.nParticlesS;     
        if(isfield(optsPhys,'HSBulk'))
            HS_f       = str2func(optsPhys.HSBulk);  
        elseif(isfield(optsNum,'HSBulk'))
            HS_f       = str2func(optsNum.HSBulk); 
        else
            HS_f       = @ZeroMap;
        end     
        eta = optsPhys.eta;
    end
    function [muSC,fnCS] = ZeroMap(h1s,h2s)
        muSC = 0;
        fnCS = 0;
    end

    function PlotRosenfeldFMT_AverageDensitiesInf(FMTShape,FMTMatrices,rho)

        figure('name','Average densities');

        nStruct = FMTMatrices.AD;    
        fields  = fieldnames(nStruct);

        noRows  = ceil(sqrt(numel(fields)));
        noCols  = ceil(numel(fields)/noRows);

        for i=1:numel(fields)
          %fields(i)
          %nStruct.(fields(i))      
            subplot(noRows,noCols,i);
            
            val = nStruct.(fields{i})*rho;
            
            plot(FMTShape.AD.Pts.y2_kv(FMTShape.AD.Pts.y1_kv == inf),val,'o'); hold on
            plot(Interp1D_AD.pts2,Interp1D_AD.InterPol*val);
            xlim([min(Interp1D_AD.pts2) max(Interp1D_AD.pts2)]);    
                        
            title(fields(i));    
        end

    %     subplot(2,2,1);
    %     FMTShape.PlotFull(FMTMatrices.AD.n2*rho,1);    
    %     title('n_2');
    %      
    % 	subplot(2,2,2);
    %     FMTShape.PlotFull(FMTMatrices.AD.n1*rho,1);
    %     title('n_1');
    %     
    %     subplot(2,2,3);
    %     FMTShape.PlotFull(FMTMatrices.AD.n1_v_1*rho,1);    
    %     title('n1_v_1');    
    %     
    %     subplot(2,2,4);
    %     FMTShape.PlotFull(FMTMatrices.AD.n1_v_2*rho,1);       
    %     title('n1_v_2');                

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
        
        y0 = -kBT*log(normalization) -VAdd;
        
        y0 = [muInit; y0];
      else
          y0 = zeros(1+N1*N2,nSpecies);
      end
    
    end
    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        mu_s  = GetExcessChemPotential(x,t,mu); 
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s);                                
    end


end