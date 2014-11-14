function data = GroupMeetingJuly_1(optsPhys,optsNum,optsPlot,name)
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
        [optsNum,optsPhys,optsPlot] = Test_GroupMeetingJuly_1();
        %TestFMT_DDFT_DiffusionHalfSpace_NSpecies();
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
    R         = diag(optsPhys.sigmaS)/2;
    fBulk     = str2func(['FexBulk_',optsNum.FexNum.Fex]);
    
    N1        = PhysArea.N(1);   N2 = PhysArea.N(2);
        
    getFex = str2func(['Fex_',optsNum.FexNum.Fex]);
    
    if(nargin<3)
        optsPlot.lineColourDDFT={'r','b','g','k','m'};
        optsPlot.doDDFTPlots=true;
    end
    
    
    HS                        = HalfSpace_FMT(PhysArea,R);
    [Pts,Diff,Int,Ind,Interp] = HS.ComputeAll(optsNum.PlotArea);
    Interp1D                  = HS.ComputeInterpolationMatrix(1,(-1:0.01:0.7)',true);
    
    %HS.AD.PlotGrid();
      
    subPts.y2_kv              = (0.:0.01:3.5)';
    subPts.y1_kv              = inf*ones(size(subPts.y2_kv));           
    IP                        = HS.AD.SubShapePts(subPts);
    Interp1D_AD.InterPol      = IP(:,HS.AD.Pts.y1_kv == inf);
    Interp1D_AD.pts1          = subPts.y1_kv; 
    Interp1D_AD.pts2          = subPts.y2_kv;
    
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
    params.Pts     = HS.Pts;     
    params.Polar   = 'cart';      
    func           = str2func(['FexMatrices_',optsNum.FexNum.Fex]);       
    IntMatrFex_2D  = DataStorage('FexData',func,params,HS); %true
    IntMatrFex_1D  = Get1DMatrices(IntMatrFex_2D,HS);
    
    IntMatrFex     = IntMatrFex_1D;
    
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
%      plotTimesV = 0:(5/optsNum.TN):5;
%      for k = 1:length(plotTimesV)
%          V =  getVAdd(y1S,y2S,plotTimesV(k),optsPhys.V1);%Vext +         
%          HS.plot(V,true,false);
%          title(['t = ',num2str(plotTimesV(k))]);
%          xlim([-5 5]);
%          ylim([0.5 5]);
%          zlim([-0.2 0]);
%          pause(0.02);        
%      end
    
  
    %****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %****************************************************************        
    tic
    
    tEq       = 0;
    markComp  = (Pts.y1_kv==inf);    
    y0        = zeros(N1*N2,nSpecies);%getInitialGuess(VAdd0);        
    y0        = y0(markComp);
    
    %PlotRosenfeldFMT_AverageDensitiesInf(HS,IntMatrFex(1),ones(size(y0)));

    %PlotRosenfeldFMT_AverageDensities(HS,IntMatrFex(1),ones(size(y0)));                   
    fprintf(1,'Computing initial condition ...');
    fsolveOpts=optimset('MaxFunEvals',2000000,'MaxIter',200000);    
    x_ic_1D        = fsolve(@f,y0,fsolveOpts);     
    rho_ic1D       = exp((x_ic_1D-Vext(markComp))/kBT);
    
    x_ig           = repmat(x_ic_1D,N1,1);
    rho_ig         = repmat(rho_ic1D,N1,1);
    
    figure; set(gcf, 'Position', [0 0 1500 1000]);    set(gcf,'Color','white');
    plot(Pts.y2_kv(markComp),rho_ic1D,'o','markersize',8,'markerFace','green'); hold on
    plot(Interp1D.pts2,Interp1D.InterPol*rho_ig,'linewidth',1.5);
    xlim([0.5 3.5]);    
    save2pdf('1DProfile.pdf',gcf);
            
    figure; set(gcf, 'Position', [0 0 2000 1500]);    set(gcf,'Color','white');
    PlotRosenfeldFMT_AverageDensitiesInf(HS,IntMatrFex(1),rho_ic1D);
    save2pdf('1DAverageDensities.pdf',gcf);
        
    %****************************************************************
    %********** Solve for equilibrium 3D final condition   **********
    %****************************************************************
    disp('Compute 3D final condition ...');
            
    markComp       = true(size(Pts.y1_kv));
    IntMatrFex     = IntMatrFex_2D;
   %x_ic           = zeros(size(x_ic));    
    fcParams             = v2struct(optsPhys,optsNum);
    fcParams.optsPhys.V1 = rmfield(fcParams.optsPhys.V1,'grav_end');
    fcParams.optsNum     = rmfield(fcParams.optsNum,'Tmax');
    fcParams.optsNum     = rmfield(fcParams.optsNum,'TN');    
    fcParams.tEq         = inf;
    
    fcParams.fsolveOpts = optimset('MaxFunEvals',2000000,'MaxIter',200000);    
    x_fc           = DataStorage('EquilibriumData_GPJuly',@ComputeEquilibriumGP,fcParams,x_ig);         	
    rho_fc         = exp((x_fc-Vext)/kBT);
        
    figure; set(gcf, 'Position', [0 0 1500 1000]);    set(gcf,'Color','white');
    HS.plot(rho_fc); title('Final Condition'); view([2,5,2]);	
    save2pdf('3DFinalCondition.pdf',gcf);
	%PlotRosenfeldFMT_AverageDensities(HS,IntMatrFex(1),rho_ic);    
    %****************************************************************
    %******** Solve for equilibrium 3D initial condition   **********
    %****************************************************************
    disp('Compute 3D initial condition ...');
        
    tEq            = 0;
    markComp       = true(size(Pts.y1_kv));
    IntMatrFex     = IntMatrFex_2D;
   %x_ic           = zeros(size(x_ic));    
    icParams             = v2struct(optsPhys,optsNum);
    icParams.optsPhys.V1 = rmfield(icParams.optsPhys.V1,'grav_end');
    icParams.optsNum     = rmfield(icParams.optsNum,'Tmax');
    icParams.optsNum     = rmfield(icParams.optsNum,'TN');
    icParams.tEq         = 0;
    
    icParams.fsolveOpts = optimset('MaxFunEvals',2000000,'MaxIter',200000);    
    x_ic           = DataStorage('EquilibriumData_GPJuly',@ComputeEquilibriumGP,icParams,x_ig);         	
    rho_ic         = exp((x_ic-Vext)/kBT);
        
    figure; set(gcf, 'Position', [0 0 1500 1000]);    set(gcf,'Color','white');

    HS.plot(rho_ic); title('Initial Condition'); view([2,5,2]);
    save2pdf('3DInitialCondition.pdf',gcf);
        
	PlotRosenfeldFMT_AverageDensities(HS,IntMatrFex(1),rho_ic);                
    save2pdf('3DAverageDensities.pdf',gcf);
    
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
    
    data     = v2struct(X_t,rho_t,mu,flux_t);%Pts,Diff,Int,Ind,Interp
    data.shape = HS;
    optsNum.plotTimes = plotTimes;  
    PlotDDFT(v2struct(optsPhys,optsNum,data));                 

   
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************                  
    function y = ComputeEquilibriumGP(params,y0)
        tEq = params.tEq;
        fprintf(1,'Computing equilibrium condition ...');        
        y          = fsolve(@f,y0,params.fsolveOpts);     
    end

    function y = f(x)
        %solves for T*log*rho + Vext                        
        y          = GetExcessChemPotential(x,tEq,mu);         
        y          = y(:);
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
        mu_s((Pts.y1_kv==inf) | (Pts.y1_kv==-inf) | (Pts.y2_kv==inf),:) = 0;
        h_s      = Diff.grad*x - Vext_grad;
        h_s([Pts.y2_kv==inf;Pts.y2_kv==inf],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==inf;false(N1*N2,1)],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==-inf;false(N1*N2,1)],:) = 0; %here, we have assumed that grad(mu) converges fast enough
                
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
                
        dxdt(Ind.bound,:)  = x(Ind.bound,:) - x_ic(Ind.bound,:);   
        %Boundary Conditions at infinity        
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
    function [muSC,fnCS] = ZeroMap(~,~)
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
            
          %  plot(FMTShape.AD.Pts.y2_kv(FMTShape.AD.Pts.y1_kv == inf),val,'o'); hold on
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