function sol = ThreePhaseContactLine_FMT_Dynamics(configIn,dirFolder)
%************************************************************************* 
%epsilon_w  | theta_CS_deg 
%--------------------------
%0.2        | 140
%0.6        | 67
%0.7        | 41.5
%0.75       | 21.5
%0.76       | 14
%*************************************************************************            
    
    disp('** ThreePhaseContactLine_FMT_DynamicsOverdamped **');    
    AddPaths();
    global PersonalUserOutput dirData
    %************************************************
    %***************  Initialization ****************
    %************************************************       
    PlotArea = struct('y1Min',-2,'y1Max',4,'y2Min',0.5,'y2Max',5,'N1',80,'N2',80);
    
    orgdirData = dirData;    
    if(nargin < 2) 
        ChangeDirData([dirData filesep 'FMT_ContactLine_Dynamics']);
    else
        ChangeDirData([dirData filesep dirFolder]);
    end    
    saveFigs   = true;
        
    diaryFile = [dirData filesep 'LogFile.txt'];
    diary(diaryFile); diary on;
       
    if((nargin > 0) && islogical(configIn))
        PhysArea = struct('N',[20,20],'L1',2,'L2',2,'y2wall',0.,...
                          'N2bound',24,'h',1,'L2_AD',2.,'alpha_deg',90);

        PhysArea.Conv  = struct('L',[],'L2',1.,'N',[50,50]);         

        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                           'Ncircle',1,'N1disc',35,'N2disc',34);

        optsNum = struct('PhysArea',PhysArea,...
                         'FexNum',Fex_Num,...
                         'maxComp_y2',inf,...
                         'y1Shift',0,...
                         'plotTimes_T',100,...
                         'plotTimes_Interval',0.1);

        V1 = struct('V1DV1','Vext_Cart_7',...
                            'epsilon_w',0.482218,...
                            'epsilon_w_max',0.3,....
                            'tau',20);
                        
        %V2 = struct('V2DV2','Phi2DLongRange','epsilon',1); 
        V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1); 

        optsPhys = struct('V1',V1,'V2',V2,...                      
                          'kBT',0.6,...   
                          'gamma',2,...
                          'Inertial',1,...
                          'Dmu',0.0,'nSpecies',1,...
                          'sigmaS',1);      
                      
        configuration = v2struct(optsNum,optsPhys);
        configName    = SaveConfig(configuration);
    elseif((nargin > 0) && isstruct(configIn))
        configuration = configIn;
        configName    = SaveConfig(configuration);        
    else
        if((nargin > 0) && (ischar(configIn)))% islogical(configIn))
            DataFolder = [dirData filesep 'Configurations' filesep];            
        else
            [configIn,DataFolder] = uigetfile([dirData filesep 'Configurations' filesep '*.mat'],['Select Config File']);
        end
        load([DataFolder,configIn]);
        disp(['Loading configuration ',[DataFolder,configIn],' ...']);
        configName = configIn(1:end-4);        
    end
    optsNum  = configuration.optsNum;
	optsPhys = configuration.optsPhys;            
    
    ChangeDirData([dirData filesep 'deg',num2str(optsNum.PhysArea.alpha_deg,3)]);
    %************************************************
    %***************  Use easier names **************
    %************************************************
    close all;            	
        
    PhysArea     = optsNum.PhysArea;           
    Fex_Num      = optsNum.FexNum;
    N1           = PhysArea.N(1);
    N2           = PhysArea.N(2);
    
    kBT             = optsPhys.kBT;    
	Dmu             = optsPhys.Dmu;    
    R               = optsPhys.sigmaS/2;     
                 
    theta_CS        = PhysArea.alpha_deg*pi/180; 
	optsPhys.HSBulk = (['FexBulk_',optsNum.FexNum.Fex]);      
	getFex          = str2func(['Fex_',optsNum.FexNum.Fex]);
    nSpecies        = 1;                     
    
    I     = eye(N1*N2);   
    eyes  = repmat(I,1,2);    
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    tic                
    %(1) Thermodynamic Values
    [rhoGas_sat,rhoLiq_sat,mu_sat] = BulkSatValues(optsPhys,[0.01;0.6;-2]);
    GetCriticalPoint(optsPhys);

    optsPhys.mu_sat      = mu_sat;
    optsPhys.rhoGas_sat  = rhoGas_sat;
    optsPhys.rhoLiq_sat  = rhoLiq_sat;
    
    %(2) Numerical Integration, Differentiation
    optsHS       = PhysArea;
    optsHS.alpha = theta_CS;
    HS                 = HalfSpace_FMT(optsHS,diag(optsPhys.sigmaS)/2);
    [Pts,Diff,Int,Ind] = HS.ComputeAll();
    HS.InterpolationPlotCart(PlotArea,true);
    PtsCart  = HS.GetCartPts();
    
    
    %(3) Numerical Convolution
	opts.V2                 = optsPhys.V2;    
    opts.nSpecies           = optsPhys.nSpecies;    
	opts.optsNum.PhysArea   = optsNum.PhysArea;    
    opts.Comments           = ['ConfigFile: ',configName,'.txt'];
    convStruct              = DataStorage(['HalfSpace_FMT' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,HS,false);   
    Conv                    = convStruct.Conv;            
        
     %(3.1) Test Convolution    
     %at infinity
     fMF = str2func(optsPhys.V2.V2DV2);
     [h1,h2,a] = fMF(1,optsPhys.V2);
     marky2Inf = (HS.Pts.y2_kv == inf);
     PrintErrorPos(Conv(marky2Inf,:)*ones(N1*N2,1)- 2*a,'Error for convolution at y2 = infinity',HS.Pts.y1_kv(marky2Inf));          
     
     %convolution profile
     if(strcmp(optsPhys.V2.V2DV2,'Phi2DLongRange'))
         y0R = PtsCart.y2_kv-R;
         h = optsPhys.V2.epsilon*(-pi^2/2 + ...
                                  +pi*atan(-y0R)+...
                                  -pi*(y0R)./(1+y0R.^2));
         PrintErrorPos(h-Conv*ones(N1*N2,1),'Error of Phi2DLongRange*1',PtsCart);                           
     end     
    %*******************************************         
    %(4) FMT Matrices
    
    fprintf(1,'Computing Fex matrices ...\n');   
        
    params.sigmaS   = optsPhys.sigmaS;
    params.FexNum   = optsNum.FexNum;
    params.PhysArea = optsNum.PhysArea;        
    params.Polar    = 'cart';      
    params.Comments = configName;
    func            = str2func(['FexMatrices_',optsNum.FexNum.Fex]);
    [IntMatrFex,recFex] = DataStorage(['HalfSpace_FMT' filesep func2str(func)],func,params,HS,false); %true      
    
    %Test
    if(recFex)
        CheckAverageDensities_Rosenfeld_3D(HS,IntMatrFex,true);
    else
        CheckAverageDensities_Rosenfeld_3D(HS,IntMatrFex,false);
    end
     
     %(5) External Potential           
    [Vext,Vext_grad]  = getVBackDVBack(PtsCart.y1_kv,PtsCart.y2_kv,optsPhys.V1);      
    VAdd     = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,0,optsPhys.V1);
    t_preprocess = toc;

    
    
    %(3) Test Integration
    a = 2;
    %(a) Integration in HalfSpace
    pts   = HS.GetCartPts();
	r     = sqrt(pts.y1_kv.^2+(pts.y2_kv-R).^2);
    testF = exp(-(r/a).^2);    
    PrintErrorPos(HS.Int*testF-a^2*pi/2,['Integration of exp(-(r/',num2str(a),')^2)']);    
    
    %(b) Integration in Composed Half Space
    pts   = HS.AD.GetCartPts();
	r     = sqrt(pts.y1_kv.^2+pts.y2_kv.^2);
    testF = exp(-(r/a).^2);
    Inth  = HS.AD.ComputeIntegrationVector();
    PrintErrorPos(Inth*testF-a^2*pi/2,['Integration of exp(-(r/',num2str(a),')^2)']);    
    
	%****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %******************** Compute Surface tension *******************
    %****************************************************************
    alpha_YCA = ComputeST();    
    if(isempty(optsNum.maxComp_y2))                
        if(theta_CS == pi/2)
            a = 3;
            optsNum.maxComp_y2 = R+a*abs(tan(alpha_YCA));            
        else
            a = 2;
            optsNum.maxComp_y2 = R+a*abs(tan(alpha_YCA)*tan(theta_CS)/(abs(tan(theta_CS))-1));
        end
    end

    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    tic    
    fprintf('Solving for equilibrium condition...\n');    
     
    p      = (rho1D_lg-rhoGas_sat)/(rhoLiq_sat-rhoGas_sat);    
    rho_ig = kron(p,rho1D_wl) + kron(1-p,rho1D_wg);         
    x_ig             = kBT*log(rho_ig)+Vext;
    
    opts             = PhysArea;
    opts.optsPhys    = optsPhys;    
    if(isfield(opts.optsPhys,'gamma'))
        opts.optsPhys    = rmfield(opts.optsPhys,'gamma');
    end
    if(isfield(opts.optsPhys,'Inertial'))
        opts.optsPhys    = rmfield(opts.optsPhys,'Inertial');
    end
    
    if(isfield(opts.optsPhys.V1,'epsilon_w_end'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_end');
    elseif(isfield(opts.optsPhys.V1,'epsilon_w_Amplitude'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_Amplitude');
    elseif(isfield(opts.optsPhys.V1,'epsilon_w_max'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'epsilon_w_max');                
    end
    if(isfield(opts.optsPhys.V1,'tau'))
        opts.optsPhys.V1 = rmfield(opts.optsPhys.V1,'tau');
    end
    opts.maxComp_y2  = optsNum.maxComp_y2;
    opts.Comments    = configName;
    
    mark             = (PtsCart.y2_kv <= optsNum.maxComp_y2);    
    [x_ic,recEq,paramsEq] = DataStorage('EquilibriumSolutions',...
                        @ComputeEquilibriumCondition,opts,x_ig); %true   
    x_ic(mark)  = x_ic;
    x_ic(~mark) = x_ig(~mark);                                       
    rho         = exp((x_ic-Vext)/kBT);
    t_eqSol     = toc;    

    %***********************************************************
    %****************  Postprocess Equilibrium  ****************
    %***********************************************************

	%Compute contact angle from density profile
    rhoV = (rhoLiq_sat+rhoGas_sat)/2; x2 = 0;
    [alphaM,pt1,pt2] = MeasureContactAngle();    
    
    if(PersonalUserOutput)        
        
        figure('Color','white','Position',[0 0 1200 800]);
        optDetails.clabel = true;        
        optDetails.nContours = [0.1,0.2,0.3,0.4,0.5,0.6,0.7];        
        HS.doPlots(rho,'contour',optDetails);  hold on;  
        plot([pt1.y1_kv;pt2.y1_kv],[pt1.y2_kv;pt2.y2_kv],'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');        
        plot((y2I-b)/slope,y2I,'--m','linewidth',1.5);
        
        if(saveFigs)       
            print2eps([dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename '_contour'],gcf);
            saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename '_contour.fig']);
        end
        
        %***************************************************************
        figure('Color','white','Position',[0 0 1200 800]);
        HS.doPlots(rho,'SC');
        zlabel('$\varrho$','Interpreter','Latex','fontsize',26);
        colormap(hsv);
        set(gca, 'CLim', [0, 1.0]);
        view([-10 5 3]);
        
        if(saveFigs)       
            print2eps([dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename],gcf);
            saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep paramsEq.Filename '.fig']);
        end        
    end
    
    if(~isfield(optsNum,'plotTimes_T'))
        sol.Measured_StaticCA = alphaM;
        
        diary    
        dirData = orgdirData;    
        return;
    end
    
    %***************************************************************    
    %***************************************************************
    %   Dynamics:
    %***************************************************************
    %***************************************************************
    
    disp('Compute dynamic 2D evolution ..');         
    mark             = true(N1*N2,1);   
    optsNum.maxComp_y2 = inf;
    opts             = PhysArea;
    opts.optsPhys    = optsPhys;     
    opts.optsNum     = optsNum;     
    opts.optsNum     = rmfield(opts.optsNum,'PhysArea');     
    opts.Comments    = configName;
    [sol,rec,paramsDyn]   = DataStorage('DynamicSolutions',@ComputeDynamics,opts,x_ic);                     
    plotTimes         = sol.outTimes;                

    %***************************************************************
    %*******************Postprocess Dynamics************************
    %***************************************************************
    
    XRho_t(mark,:)  = (sol.X_t(:,1:N1*N2))';
	rho_t           = exp((XRho_t-Vext*ones(1,length(plotTimes)))/kBT);
    if(optsPhys.Inertial)         
         flux_t      = [rho_t;rho_t].*(sol.X_t(:,1+N1*N2:end))';
    else
        flux_t          = zeros(2*N1*N2,length(plotTimes));    
        for i = 1:length(plotTimes)
            flux_t(:,i) = GetFlux(XRho_t(:,i),plotTimes(i));
        end
    end            

    shapeSub      = struct('y1Min',-2,'y1Max',2,'y2Min',R,'y2Max',3,'N',[40,40]);
    subBox        = Box(shapeSub);
    IP            = HS.SubShapePtsCart(subBox.Pts);
    Int_SubOnFull = subBox.ComputeIntegrationVector()*IP;
        
    Int_of_path   = subBox.IntFluxThroughDomain(100)*blkdiag(IP,IP);
    
    accFlux = zeros(length(plotTimes),1);    
    for i = 2:length(plotTimes)        
        accFlux(i) = accFlux(i-1)+0.5*(plotTimes(i)-plotTimes(i-1))*...
            (Int_of_path*(flux_t(:,i-1)+flux_t(:,i)));
    end    
    
    Subspace      = v2struct(Int_of_path,Int_SubOnFull,accFlux);%InterpPath,Path
    Subspace.subArea = subBox;
    %***************************************************************
    %***********************Plotting********************************
    %***************************************************************
    optsPlot.nContours = [0.1,0.3,0.5,0.7];                
    data       = v2struct(rho_t,flux_t,optsPlot,Subspace);
    data.shape = HS;
    optsNum.plotTimes = plotTimes; 

    plotData = v2struct(optsPhys,optsNum,data);
         
    if(PersonalUserOutput)        
        PlotDDFT(plotData,[dirData filesep 'DynamicSolutions' filesep paramsDyn.Filename(1:end-4)]);            
        PlotDDFT(plotData);
    end
    
    diary    
    dirData = orgdirData;
    
    %***************************************************************
    %***********************End of main code************************
    %***************************************************************
    
        
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************      
    function sol = ComputeDynamics(params,x_ic)
        plotTimes  = (0:optsNum.plotTimes_Interval:optsNum.plotTimes_T);
        if(optsPhys.Inertial)
            mM            = ones(3*N1*N2,1);    
            mM([Ind.top;false(N1*N2,1);Ind.bottom]) = 0;            
            mM            = mM([mark;mark;mark]);
            opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
            [sol.outTimes,sol.X_t] =  ode15s(@dx_dtInertial,plotTimes,...
                                [x_ic(mark);zeros(2*sum(mark),1)],opts);
        else
            mM             = ones(N1*N2,1);    
            mM(Ind.bottom|Ind.top) = 0;
            mM             = mM(mark);
            opts = odeset('RelTol',10^-10,'AbsTol',10^-10,'Mass',diag(mM));
            [sol.outTimes,sol.X_t] =  ode15s(@dx_dt,plotTimes,x_ic(mark),opts);                
        end
    end    
    function [x_ic] = ComputeEquilibriumCondition(params,x_ig)               
        [x_ic,~,flag_fsolve]    = fsolve(@GetExcessChemPotential,x_ig(mark));
        if(flag_fsolve~=1)
            x_ic = 0;
        end
    end
    function mu_s = GetExcessChemPotential(xm)        
        x(mark)  = xm;
        x(~mark) = x_ig(~mark);
        x        = x';
        mu_s = GetExcessChemPotentialFull(x);
        mu_s = mu_s(mark);
    end   
    function mu_s = GetExcessChemPotentialFull(x)                
        rho_s    = exp((x-Vext)/kBT);
        mu_s     = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - (mu_sat + Dmu);
        end
        mu_s = mu_s + x + Conv*rho_s + VAdd;           
    end   
    function dxdt = dx_dt(t,x)
       
        VAdd     = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,t,optsPhys.V1);                
        
        x(mark)  = x;
        x(~mark) = x_ic(~mark);
        %x = x';
        
        mu_s     = GetExcessChemPotentialFull(x);        
        mu_s((Pts.y2_kv==inf),:) = 0;
        h_s      = Diff.grad*x - Vext_grad;
        h_s([Pts.y2_kv==inf;Pts.y2_kv==inf],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==inf;false(N1*N2,1)],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s([Pts.y1_kv==-inf;false(N1*N2,1)],:) = 0; %here, we have assumed that grad(mu) converges fast enough
                
        dxdt     = kBT*Diff.Lap*mu_s + eyes*(h_s.*(Diff.grad*mu_s));  
        dxdt     = dxdt/optsPhys.gamma;
                
        %Boundary Conditions at infinity
        flux_dir           = Diff.grad*mu_s;
        dxdt(Ind.bottom,:) = Ind.normalBottom*flux_dir;
        dxdt(Ind.top,:)    = x(Ind.top)-x_ic(Ind.top);

        dxdt = dxdt(:);      
        dxdt = dxdt(mark);
    end
    function dxdt = dx_dtInertial(t,x)         
        VAdd     = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,t,optsPhys.V1);                
        
        x        = reshape(x,sum(mark),3);
        x(mark,:)  = x;
        x(~mark,1) = x_ic(~mark);  %outside computational domain, set density to intial and set velocity to zero
        
        y        = x(:,1);
        u        = x(:,2);
        v        = x(:,3);
        uv       = x(:,2:3);
        uv       = uv(:);
        rho      = exp((y-Vext)/kBT);        
        
        %Convective term:        
        %C        = [diag(u) diag(v)]*Diff.grad;
        %C        = sparse(diag(u))*Diff.Dy1 + sparse(diag(v))*Diff.Dy2;
        m        = N1*N2;
        C        = sparse(1:m,1:m,u)*Diff.Dy1 + sparse(1:m,1:m,v)*Diff.Dy2;
        %CC      = blkdiag(C,C);
        mu_s     = GetExcessChemPotentialFull(y);
        mu_s(~mark) = 0;
        mu_s((Pts.y2_kv==inf),:) = 0; %necessary?
        
        %Equations       
        h_s1    = Diff.Dy1*y - Vext_grad(1:m);
        h_s2    = Diff.Dy2*y - Vext_grad(1+m:end);
                
        h_s1(Pts.y2_kv==inf)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s2(Pts.y2_kv==inf)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s2(Pts.y1_kv==inf)  = 0; %here, we have assumed that grad(mu) converges fast enough
        h_s2(Pts.y1_kv==-inf) = 0; %here, we have assumed that grad(mu) converges fast enough
        
        %h_s      = Diff.grad*y - Vext_grad;
        %h_s([Pts.y2_kv==inf;Pts.y2_kv==inf],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        %h_s([Pts.y1_kv==inf;false(N1*N2,1)],:)  = 0; %here, we have assumed that grad(mu) converges fast enough
        %h_s([Pts.y1_kv==-inf;false(N1*N2,1)],:) = 0; %here, we have assumed that grad(mu) converges fast enough
        
        dydt     = - kBT*(Diff.div*uv) - (u.*h_s1  + v.*h_s2);        
        %dydt     = - kBT*(Diff.div*uv) - [diag(u) diag(v)]*h_s;        
        %duvdt    = - CC*uv - optsPhys.gamma*uv - Diff.grad*mu_s;
        duvdt    = - [C*u;C*v] - optsPhys.gamma*uv - Diff.grad*mu_s;
               
        %Boundary Conditions for velocities
        dydt(Ind.top)                       = y(Ind.top)- x_ic(Ind.top); %Density at y_2 = infinity        
        duvdt([false(N1*N2,1);Ind.bottom])  = uv([false(N1*N2,1);Ind.bottom]); %Normal Velocity at wall       
        dxdt = [dydt(mark); duvdt([mark;mark])];       
    end
    function flux = GetFlux(x,t)
        rho_s = exp((x-Vext)/kBT);       
        VAdd  = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,t,optsPhys.V1);                
        mu_s  = GetExcessChemPotentialFull(x); 
        mu_s((Pts.y2_kv==inf),:) = 0;
        flux  = -[rho_s;rho_s].*(Diff.grad*mu_s)/optsPhys.gamma;                                
    end
    function z = rhoX1(x1)
        IP = HS.ComputeInterpolationMatrix(x1,x2);
        z  = IP.InterPol*rho-rhoV;
    end
    function alpha_YCA = ComputeST()
        optss            = optsPhys;                
        optss.rho_iguess = (rhoLiq_sat+rhoGas_sat)/2 + ...
                              (rhoLiq_sat-rhoGas_sat)/2*tanh((Pts.y1-optsNum.y1Shift)*sin(theta_CS));
                          
        [rho1D_lg,parms] = FMT_1D_Interface(HS,IntMatrFex,optss,Fex_Num,Conv,false,optsNum.y1Shift);
        om_LiqGas        = parms.Fex;

        optss.rho_iguess = rhoGas_sat;
        [rho1D_wg,parms] = FMT_1D(HS,IntMatrFex,optss,Fex_Num,Conv,false);
        om_wallGas       = parms.Fex;

        optss.rho_iguess = rhoLiq_sat;
        [rho1D_wl,parms] = FMT_1D(HS,IntMatrFex,optss,Fex_Num,Conv,false);
        om_wallLiq       = parms.Fex;
        
        close all;    

        
        fprintf(['Omega(Liq/Gas) = ',num2str(om_LiqGas),'\n']);
        fprintf(['Omega(wall/Liq) = ',num2str(om_wallLiq),'\n']);
        fprintf(['Omega(wall/Gas) = ',num2str(om_wallGas),'\n']);
        

        alpha_YCA = ComputeContactAngle(om_wallGas,om_wallLiq,om_LiqGas);
    end

    function  [alphaM,pt1,pt2] = MeasureContactAngle()
    
        upperY2_angleComp = min(optsNum.maxComp_y2,15);        
        fsolveOpts=optimset('Display','off');    
        x2_1 = HS.CompSpace2(upperY2_angleComp/(2*sin(theta_CS)));
        x2   = x2_1;
        x1_1 = fsolve(@rhoX1,-0.5,fsolveOpts);
        x2_2 = HS.CompSpace2(upperY2_angleComp/sin(theta_CS));
        x2   = x2_2;
        x1_2 = fsolve(@rhoX1,-0.5,fsolveOpts);

        [y1_1,y2_1] = HS.PhysSpace(x1_1,x2_1);
        [y1_2,y2_2] = HS.PhysSpace(x1_2,x2_2);
        pt1 = HS.GetCartPts(y1_1,y2_1);
        pt2 = HS.GetCartPts(y1_2,y2_2);

        slope  = (pt2.y2_kv-pt1.y2_kv)/(pt2.y1_kv-pt1.y1_kv);
        b      = pt2.y2_kv - slope*pt2.y1_kv;
        alphaM = mod(atan(slope),pi);
        y2I    = (PlotArea.y2Min:0.1:PlotArea.y2Max)';
        fprintf(['Measured Contact Angle: ',num2str(alphaM*180/pi),' [deg]\n']);
    end
end