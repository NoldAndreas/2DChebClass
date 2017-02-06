function ThreePhaseContactLine_FMT_HalfSpaceSkewed(CompCase,Case90)
%************************************************************************* 
% solve
% 0 = kBT*log(rho) + int(rho(r')*Phi2D(r-r'),dr') +...
%               + mu_HS(rho_s) + V_ext - mu_sat
%
%************************************
%epsilon_w  | theta_CS_deg 
%--------------------------
%0.2        | 140
%0.6        | 67
%0.7        | 41.5
%0.75       | 21.5
%0.76       | 14
%*************************************************************************            
    
    disp('** ThreePhaseContactLine FMT HalfSpaceSkewed **');    
    AddPaths();
    global dirData
    %************************************************
    %***************  Initialization ****************
    %************************************************        
    if(nargin == 0)
        CompCase = 4;
    end        
    
    orgdirData = dirData;
    
    dirData = [dirData filesep 'FMT_ContactLine_Version1_2013Oct12' filesep];
    saveFigs = true;
    
    if((nargin == 2) && Case90)
        theta_CS_deg = 90;
        dirData      = [orgdirData filesep 'FMT_ContactLine_Version1_2013Oct12' filesep 'deg90'];            
    else
        Case90 = false;
    end
    
    plotX0 = 2;
       
    switch CompCase
        case 1
            epsilon_w    = 0.76;
            if(Case90)
                y2MaxComp    = 2;
                plotX0 = 5;
            else
                theta_CS_deg = 14;
                y2MaxComp    = 4;            
                dirData      = [dirData 'deg14'];
            end
        case 2
            epsilon_w    = 0.75;
            if(Case90)
                y2MaxComp    = 2;
                plotX0 = 5;
            else
                theta_CS_deg = 21.5;
                y2MaxComp    = 8;            
                dirData      = [dirData 'deg21_5'];
            end            
        case 3
            epsilon_w    = 0.7;            
            if(Case90)
                y2MaxComp    = 3.5;
            else
                theta_CS_deg = 41.5;                
                y2MaxComp    = 20;
                dirData      = [dirData 'deg41_5'];
            end
        case 4
            epsilon_w    = 0.6;            
            if(Case90)
                y2MaxComp    = 5;
            else
                theta_CS_deg = 67;               
                y2MaxComp    = 20;
                dirData      = [dirData 'deg67'];
            end
        case 5
            epsilon_w    = 0.2;
            if(Case90)
                y2MaxComp    = 5;
            else
                theta_CS_deg = 140;                
                y2MaxComp    = 20;
                dirData      = [dirData 'deg140'];
            end
        otherwise
            epsilon_w     = 0.5;            
            theta_CS_deg  = 90;            
            y2MaxComp     = 20;
            dirData       = [dirData 'deg90'];            
    end
    
    y2MaxCompPlot = 6;

    if(~exist(dirData,'dir'))
        mkdir(dirData);
    end
    
    diaryFile = [dirData filesep 'LogFile.txt'];
    if(~exist(diaryFile,'file'))
        fid = fopen(diaryFile,'w');
        fclose(fid);
    end
    diary(diaryFile); diary on;
       
    N1           = 40;
    N2           = 40;    
    theta_CS     = theta_CS_deg*pi/180; 
        
	PhysArea = struct('N',[N1,N2],'L1',2/sin(theta_CS),'L2',2,'y2wall',0.,...
                      'N2bound',24,'h',1,'L2_AD',2.,'alpha_deg',theta_CS_deg);
                   
    PhysArea.Conv  = struct('L',[],'L2',1.,'N',[50,50]);
    
    Plot_Area = struct('y1Min',min(-plotX0,-plotX0+y2MaxCompPlot/tan(theta_CS)),...
                       'y1Max',max(plotX0,plotX0+y2MaxCompPlot/tan(theta_CS)),...
                       'N1',80,'N2',80,...
                       'y2Min',0.5,'y2Max',y2MaxCompPlot+2);
    
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',35,'N2disc',34);
        
    optsNum = struct('PhysArea',PhysArea,...
                     'PlotArea',Plot_Area,...                     
                     'FexNum',Fex_Num,...
                     'maxComp_y2',y2MaxComp,...
                     'y1Shift',0,...%-0.5
                     'DDFTCode','DDFT_HalfSpace_FMT_2Phase_Sat');
                      
    V1 = struct('V1DV1','Vext_Cart_7','epsilon_w',epsilon_w);%0.4);
    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1); 
                     
    optsPhys = struct('V1',V1,'V2',V2,...                      
                      'kBT',0.6,...   
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);        
    %************************************************
    %***************  Initialization ****************
    %************************************************
    close all;
        
    PhysArea     = optsNum.PhysArea;           
    kBT          = optsPhys.kBT;    
    R            = optsPhys.sigmaS/2;     
    
	optsPhys.HSBulk = (['FexBulk_',optsNum.FexNum.Fex]);      
	getFex          = str2func(['Fex_',optsNum.FexNum.Fex]);
    nSpecies        = 1;
            
    Phi_r          = str2func(optsPhys.V2.V2DV2);
    Dmu            = optsPhys.Dmu;
    
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
    HS.InterpolationPlotCart(optsNum.PlotArea,true);
    PtsCart  = HS.GetCartPts();
    
    
    %(3) Numerical Convolution
	opts.V2                 = optsPhys.V2;    
    opts.nSpecies           = optsPhys.nSpecies;    
	opts.optsNum.PhysArea   = optsNum.PhysArea;    
    [convStruct,recomputed] = DataStorage(['HalfSpace_FMT' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,HS);   
    Conv                    = convStruct.Conv;            
        
     %(3.1) Test Convolution
    
     %(a) Test Interpolation for Convolution with current function
    % yPtsCheck          = [0 2 ; 0 PhysArea.y2wall ; -10 0 ; 20 10;0 Inf];
	% if(recomputed)        
    %    HS.TestConvolutionMatrix(yPtsCheck,@Phi);
    % else
    %     HS.TestConvolutionMatrix(yPtsCheck,@Phi,false);
    % end
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
         
         
         %test convolution with a function at infinity         
         testF = 1./(1+(Pts.y1_kv));
     end
     
     %(b) Test complute convolution computation with toy function 1/((1+x^2+y^2)^2)
%     opts.PhysArea = optsNum.PhysArea;
%     opts.FexNum   = optsNum.FexNum;
%     ConvTest      = DataStorage('ComputeConvolutionMatrix',@ComputeTestConvMatrix,opts,HS);          
%     h             = convTestResult(PtsCart.y1_kv,PtsCart.y2_kv);
%     PrintErrorPos(h-ConvTest*ones(N1*N2,1),'Test Convolution "{1/((1+x^2+y^2)^2)}*1"',HS.Pts);          
     %*******************************************    
     
      %(4) FMT Matrices
    
    fprintf(1,'Computing Fex matrices ...\n');   
        
    params.sigmaS  = optsPhys.sigmaS;
    params.FexNum  = optsNum.FexNum;
    params.PhysArea = optsNum.PhysArea;
        
    params.Polar   = 'cart';      
    func           = str2func(['FexMatrices_',optsNum.FexNum.Fex]);                
    [IntMatrFex,recFex] = DataStorage(['HalfSpace_FMT' filesep func2str(func)],func,params,HS); %true      
     
     %(5) External Potential           
    Vext     = getVBackDVBack(PtsCart.y1_kv,PtsCart.y2_kv,optsPhys.V1);      
    VAdd     = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,0,optsPhys.V1);
    t_preprocess = toc;

    %(2) Test FMT
     if(recFex)
        CheckAverageDensities_Rosenfeld_3D(HS,IntMatrFex,true);
    else
        CheckAverageDensities_Rosenfeld_3D(HS,IntMatrFex,false);
    end
    
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
    
    fprintf('*************************\n');
    %********************************
    
    
	%****************************************************************
    %**************** Solve for equilibrium 1D condition   **********
    %****************************************************************        
    
    %***** Compute Surface tension *****    
    optss            = optsPhys;        
    optss.Dmu        = optsPhys.Dmu; 

    olddirData  = dirData;
    dirData     = [dirData filesep 'Case_' num2str(CompCase)];

    optss.rho_iguess = (rhoLiq_sat+rhoGas_sat)/2 + ...
                          (rhoLiq_sat-rhoGas_sat)/2*tanh((Pts.y1-optsNum.y1Shift)*sin(theta_CS));
    [rho1D_lg,parms] = FMT_1D_Interface(HS,IntMatrFex,optss,Fex_Num,Conv,true,optsNum.y1Shift);
    om_LiqGas        = parms.Fex;

    optss.rho_iguess = rhoGas_sat;
    [rho1D_wg,parms] = FMT_1D(HS,IntMatrFex,optss,Fex_Num,Conv,true);
    om_wallGas       = parms.Fex;

    optss.rho_iguess = rhoLiq_sat;
    [rho1D_wl,parms] = FMT_1D(HS,IntMatrFex,optss,Fex_Num,Conv,true);
    om_wallLiq       = parms.Fex;

    fprintf(['Omega(Liq/Gas) = ',num2str(om_LiqGas),'\n']);
    fprintf(['Omega(wall/Liq) = ',num2str(om_wallLiq),'\n']);
    fprintf(['Omega(wall/Gas) = ',num2str(om_wallGas),'\n']);

    ComputeContactAngle(om_wallGas,om_wallLiq,om_LiqGas);
    dirData = olddirData;
    
    close all;    
    
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
    mark = (PtsCart.y2_kv < optsNum.maxComp_y2);
    opts.maxComp_y2  = optsNum.maxComp_y2;
        
    x_ic        = DataStorage(['ThreePhaseContactLine_FMT_HalfSpaceSkewed' filesep 'EquilibriumSolutions'],...
                        @ComputeEquilibriumCondition,opts,x_ig); %true   
    close all;
    x_ic(mark)  = x_ic;
    x_ic(~mark) = x_ig(~mark);
    rho         = exp((x_ic-Vext)/kBT);
    t_eqSol   = toc;    

    %************************************************
    %****************  Postprocess  ****************
    %************************************************                
    
	%***************************************************************
	%Compute contact angle from density profile
    rhoV = (rhoLiq_sat+rhoGas_sat)/2;
    fsolveOpts=optimset('Display','off');    
    x2_1 = HS.CompSpace2(optsNum.maxComp_y2/(2*sin(theta_CS)));
    x2   = x2_1;
    x1_1 = fsolve(@rhoX1,-0.5,fsolveOpts);
    x2_2 = HS.CompSpace2(optsNum.maxComp_y2/sin(theta_CS));
    x2   = x2_2;
    x1_2 = fsolve(@rhoX1,-0.5,fsolveOpts);

    [y1_1,y2_1] = HS.PhysSpace(x1_1,x2_1);
    [y1_2,y2_2] = HS.PhysSpace(x1_2,x2_2);
    pt1 = HS.GetCartPts(y1_1,y2_1);
    pt2 = HS.GetCartPts(y1_2,y2_2);

    slope  = (pt2.y2_kv-pt1.y2_kv)/(pt2.y1_kv-pt1.y1_kv);
    b      = pt2.y2_kv - slope*pt2.y1_kv;
    alphaM = mod(atan(slope),pi);
    y2I = (Plot_Area.y2Min:0.1:Plot_Area.y2Max)';
    fprintf(['Measured Contact Angle: ',num2str(alphaM*180/pi),' [deg]\n']);
    %***************************************************************
    
    dirData     = [dirData filesep 'Case_' num2str(CompCase)];

    figure('Color','white','Position',[0 0 1200 800]);
    %optDetails.nContours = 10;       
    optDetails.clabel = true;        
    optDetails.nContours = [0.1,0.2,0.3,0.4,0.5,0.6,0.7];        
    %HS.PlotGridLines(); hold on;
    HS.plot(rho,'contour',optDetails);

    hold on;  
    plot([pt1.y1_kv;pt2.y1_kv],[pt1.y2_kv;pt2.y2_kv],'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');        
    plot((y2I-b)/slope,y2I,'--m','linewidth',1.5);

    if(saveFigs)       
        print2eps([dirData filesep 'DensityPlotCL'],gcf);
        saveas(gcf,[dirData filesep 'DensityPlotCL.fig']);
    end

    %***************************************************************
    figure('Color','white','Position',[0 0 1200 800]);
    y2MaxCompPlot = 2;
    PlotArea2 = struct('y1Min',min(-plotX0,-plotX0+y2MaxCompPlot/tan(theta_CS)),...
                       'y1Max',max(plotX0,plotX0+y2MaxCompPlot/tan(theta_CS)),...
                       'N1',60,'N2',60,...                           
                       'y2Min',0.5,'y2Max',y2MaxCompPlot+2);
    HS.InterpolationPlotCart(PlotArea2,true);        
    HS.plot(rho,'SC');
    zlabel('$\varrho$','Interpreter','Latex','fontsize',26);
    colormap(hsv);
    set(gca, 'CLim', [0, 1.0]);
    view([-10 5 3]);                       

    if(saveFigs)       
        print2eps([dirData filesep 'ContourPlot'],gcf);
        saveas(gcf,[dirData filesep 'Density_Wall_AverageDensities.fig']);
    end        

    
    display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess),'\n']);
    display(['Equilibrium, Computation time (sec): ', num2str(t_eqSol),'\n']);
    
    diary
    
    dirData = orgdirData;
        
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************  
    
    function x_ic = ComputeEquilibriumCondition(params,x_ig)        
       % options = optimset('MaxFunEvals',300000);
        x_ic    = fsolve(@GetExcessChemPotential,x_ig(mark));%,options);
    end


    function mu_s = GetExcessChemPotential(xm)
        
        x(mark)  = xm;
        x(~mark) = x_ig(~mark);
        x        = x';
                
        rho_s    = exp((x-Vext)/kBT);
        mu_s     = getFex(rho_s,IntMatrFex,kBT,R);
                       
        for iSpecies=1:nSpecies
           mu_s(:,iSpecies) = mu_s(:,iSpecies) - (mu_sat + Dmu);
        end
        mu_s = mu_s + x + Conv*rho_s + VAdd;   
        mu_s = mu_s(mark);
    end   

    function z = Phi(r)
         z = Phi_r(r,optsPhys.V2);
    end       

    function z = f_ConvTest(r)
        z = 1./((1+r.^2).^2);
    end

    function z = convTestResult(y1,y2)
        y2 = y2-R;
        z  = pi/2*(1 + y2./sqrt(1+y2.^2));
        z(y2==inf) = pi;
    end
    function ConvTest = ComputeTestConvMatrix(h1,h2)
        ConvTest = HS.ComputeConvolutionMatrix(@f_ConvTest,optsNum.FexNum);
    end

    function z = rhoX1(x1)
        IP = HS.ComputeInterpolationMatrix(x1,x2);
        z  = IP.InterPol*rho-rhoV;
    end

end