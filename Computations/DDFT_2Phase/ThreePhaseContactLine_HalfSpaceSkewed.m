function ThreePhaseContactLine_HalfSpaceSkewed()
%************************************************************************* 
% data = ThreePhaseContactLine_HalfSpaceSkewed
%
% Set mu_s :=  = kBT*log(rho) + int(rho(r')*Phi2D(r-r'),dr') +...
%               + mu_HS(rho_s) + V_ext - mu_sat
%
% Equilibrium:
% (EQ 1) mu_s = 0  
% Dynamics: 
% (DYN 1) drho/dt = div(rho*grad(mu_s))
% (BC )    0      = n*grad(mu_s) 
% where n is the normal to the wall.
%*************************************************************************            

%epsilon_w   |  contact angle [deg]
%0.4         |  63
%0.5         |  33.1
%0.52        |  26.5
    %************************************************
    %*************** Parameters *********************
    %************************************************
    disp('** ThreePhaseContactLine_HalfSpaceSkewed **');
    
    close all;
    
    N1 = 40;    
    N2 = 40;
    theta_CS_A = 26.5;%63.;
    theta_CS = theta_CS_A*pi/180;
  %  theta_CS = pi/2;

    PhysArea  = struct('N',[N1,N2],'y2Min',0,'L1',4/sin(theta_CS),'L2',2,'alpha',theta_CS);%sin(theta_CS)

    Plot_Area = struct('y1Min',-10/sin(theta_CS),'y1Max',10/sin(theta_CS),'N1',100,'N2',100,...
                       'y2Min',0,'y2Max',20);

    Sub_Area = struct('y1Min',-1,'y1Max',1,'N',[20,20],...
                      'y2Min',0,'y2Max',2);    

    PhysArea.Conv  = struct('L',2,'L2',1.,'N',[50,50]);

    optsNum = struct('PhysArea',PhysArea,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.1:6,...                     
                     'DDFTCode','DDFT_DiffusionHalfSpace_2Phase_Sat');                 

    V1 = struct('V1DV1','Vext_Cart_7','epsilon_w',0.52);
    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);

    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','MuCarnahanStarling',...
                      'kBT',0.7,...                      
                      'Dmu',0.,'nSpecies',1);   
 
    %************************************************
    %***************  Initialization ****************
    %************************************************    
        
    PhysArea     = optsNum.PhysArea;   
    kBT          = optsPhys.kBT;     
        
    HS_f         = str2func(optsPhys.HSBulk);          
    Phi_r        = str2func(optsPhys.V2.V2DV2);
    Dmu          = optsPhys.Dmu;    
    	
    %************************************************
    %****************  Preprocess  ****************
    %************************************************       
    tic        
    %(1) Thermodynamic Values
    [rhoGas_sat,rhoLiq_sat,mu_sat] = BulkSatValues(optsPhys,[0.01;0.6;-2]);
    mu = mu_sat + Dmu;
    
    %(2) Numerical Integration, Differentiation
	HS                 = HalfSpaceSkewed(PhysArea);
    [Pts,Diff,Int,Ind] = HS.ComputeAll(optsNum.PlotArea);  
    
    %(3) Numerical Convolution
	opts.V2            = optsPhys.V2;    
    opts.nSpecies      = optsPhys.nSpecies;      
	opts.optsNum.PhysArea  = optsNum.PhysArea;    
    convStruct         = DataStorage(['HalfSpace' filesep 'FexMatrices_Meanfield'],@FexMatrices_Meanfield,opts,HS);   
    Conv               = convStruct.Conv;    
    %yPtsCheck          = [20 10 ; 0 2 ; 0 PhysArea.y2Min ; -10 0 ];
    %HS.TestConvolutionMatrix(yPtsCheck,@Phi);
    
    %(4) External Potential
    PtsCart      = HS.GetCartPts();
    Vext         = getVBackDVBack(PtsCart.y1_kv,PtsCart.y2_kv,optsPhys.V1);              
    Vadd         = getVAdd(PtsCart.y1_kv,PtsCart.y2_kv,0,optsPhys.V1);
    t_preprocess = toc;
    
    %****************************************************************
    %**************** Solve for 1D equilibrium condition   **********
    %****************************************************************    
    markInf = (Pts.y1_kv == inf);
    Conv1D  = Conv(markInf,markInf);
    Vext1D  = Vext(markInf);
    Vadd1D  = Vadd(markInf);
    y2_1D   = PtsCart.y2_kv(markInf);    
    x2Interp1D = (-1:0.01:HS.CompSpace2(10))';
    Interp1D = HS.ComputeInterpolationMatrix(1,x2Interp1D,true);
    
    [h1,h2,Int_1D]  = HS.ComputeIntegrationVector();
    %the values we are trying to compute here decay algebraically. 
    %Ignore values for y2>50!
    Int_1D(y2_1D > 50) = 0;
    
    y2Interp1D = sin(theta_CS)*Interp1D.pts2;
    Interp1D   = Interp1D.InterPol(:,markInf);
    
    
    %**********************************************
    %(1) Compute Wall-Gas Interface
    x1D      = fsolve(@GetExcessChemPotential1D,kBT*log(rhoGas_sat*ones(N2,1)));
    rho1D_wg = exp((x1D-Vext1D)/kBT);
    [om_ex_wg,om_ex_wgLoc] = GetExcessGrandPotential(rho1D_wg);
    disp(Int_1D*(rho1D_wg-rho1D_wg(end)));
    
    %**********************************************
    %(2) Compute Wall-Liq Interface
    x1D      = fsolve(@GetExcessChemPotential1D,kBT*log(rhoLiq_sat*ones(N2,1)));
    rho1D_wl = exp((x1D-Vext1D)/kBT);
    [om_ex_wl,om_ex_wlLoc] = GetExcessGrandPotential(rho1D_wl);
    
    %figure;
    %plot(y2_1D,om_ex_wgLoc,'o'); hold on; plot(y2Interp1D,Interp1D*om_ex_wgLoc); 
    %plot(y2_1D,-om_ex_wlLoc,'ro'); hold on; plot(y2Interp1D,-Interp1D*om_ex_wlLoc,'r'); 
    %xlabel('$y_2$','Interpreter','Latex'); ylabel('$\Omega_{ex}$','Interpreter','Latex');
    %xlim([min(y2Interp1D) max(y2Interp1D)]);
    
    figure; subplot(1,2,1);
    plot(y2_1D,rho1D_wg,'o'); hold on; plot(y2Interp1D,Interp1D*rho1D_wg); 
    plot(y2_1D,rho1D_wl,'ro'); hold on; plot(y2Interp1D,Interp1D*rho1D_wl,'r'); 
    xlabel('$y_2$','Interpreter','Latex'); ylabel('$\rho$','Interpreter','Latex');
    xlim([min(y2Interp1D) max(y2Interp1D)]);
    
    %ftest = exp(-PtsCart.y2_kv(markInf));
    %disp(['Test Integration: ',num2str(Int_1D*ftest)]);
    
    %**********************************************
    %(3) Compute Liq-Gas Interface
    markInf = (Pts.y2_kv == inf);
    Conv1D  = Conv(markInf,markInf);
    Vext1D  = Vext(markInf);
    Vadd1D  = Vadd(markInf);    
    y1_1D   = Pts.y1_kv(markInf);    
    x1Interp1D = (HS.CompSpace1(-10/sin(theta_CS)):0.01:HS.CompSpace1(10/sin(theta_CS)))';
    Interp1D = HS.ComputeInterpolationMatrix(x1Interp1D,1,true);
    [h1,Int_1D,h2] = HS.ComputeIntegrationVector(); 
    %Int_1D(abs(y1_1D)>50) = 0;
    Int_1D       = sin(theta_CS)*Int_1D; %for integration, take into consideration that integration is normal to the interface!!
    y1Interp1D   = Interp1D.pts1;
    Interp1D     = Interp1D.InterPol(:,markInf);
    
    rho_ig   = (rhoLiq_sat+rhoGas_sat)/2 + (rhoLiq_sat-rhoGas_sat)/2*tanh(y1_1D*sin(theta_CS)/3);
    
    x_ig     = kBT*log(rho_ig);
    
    x1D      = fsolve(@f_1DLG,x_ig);
    rho1D_lg = exp((x1D-Vext1D)/kBT);
    [om_ex_lg,floc] = GetExcessGrandPotential(rho1D_lg);
        
    subplot(1,2,2);
    y1_normalP = sin(theta_CS)*y1_1D;
    y1_normal  = sin(theta_CS)*y1Interp1D;
    plot(y1_normalP,rho1D_lg,'o'); hold on; plot(y1_normal,Interp1D*rho1D_lg);  
    xlim([min(y1_normal) max(y1_normal)]);    
    
    
    figure('Name','InitialGuess');
    subplot(1,2,1);
    plot(y1_normalP,rho_ig,'o'); hold on; plot(y1_normal,Interp1D*rho_ig);  
    xlim([min(y1_normal) max(y1_normal)]);
    
    subplot(1,2,2);
    plot(y1_normalP,Conv(markInf,markInf)*rho_ig,'o'); hold on; plot(y1_normal,Interp1D*(Conv(markInf,markInf)*rho_ig));
    xlim([min(y1_normal) max(y1_normal)]);
    title('Convolution of initial guess');
    
    figure('Name','Excess grand potential');
    plot(y1_normalP,Int_1D'.*floc,'o'); hold on; plot(y1_normal,Interp1D*(Int_1D'.*floc));
    xlim([min(y1_normal) max(y1_normal)]);
    title('LG Surface tension');
    
	figure('Name','Test Convolution 1');
    fT = exp(-(y1_normalP/3).^2);
    plot(y1_normalP,Conv(markInf,markInf)*fT,'o'); hold on; plot(y1_normal,Interp1D*(Conv(markInf,markInf)*fT));
    xlim([min(y1_normal) max(y1_normal)]);
    title('Convolution with Gaussion');
    
    
    herr = Int_1D*exp(-y1_normalP.^2) - sqrt(pi);
    PrintErrorPos(herr,'Integration Error 2');
    
    %Contact Angle:
    fprintf(['Omega(Liq/Gas) = ',num2str(om_ex_lg),'\n']);
    fprintf(['Omega(wall/Liq) = ',num2str(om_ex_wl),'\n']); 
    fprintf(['Omega(wall/Gas) = ',num2str(om_ex_wg),'\n']);
        
    
    theta = ComputeContactAngle(om_ex_wg,om_ex_wl,om_ex_lg);
    disp(['The contact angle is: ',num2str(theta*180/pi),'°']);
    
    %****************************************************************
    %**************** Solve for equilibrium condition   ************
    %****************************************************************    
    tic    
    fprintf('Solving for equilibrium condition...');    
    
    rho_ig = kron(rho1D_lg,ones(N2,1));
    figure('Name','Initial Guess');
    HS.doPlots(rho_ig,'SC');
    x_ig = kBT*log(rho_ig)+Vext;
    
    opts             = PhysArea;
    opts.optsPhys    = optsPhys;
    opts.ig          = x_ig;
    mark = (PtsCart.y2_kv < 20);
    opts.mark        = mark;
    
    xm_ic    = DataStorage('EquilibriumSolutions',@ComputeEquilibriumCondition,opts,HS);
	
    x_ic(mark)   = xm_ic;
	x_ic(~mark)  = x_ig(~mark);
    x_ic         = x_ic';
    
    %x_ic    = fsolve(@f,kBT*log(rhoGas_sat*ones(N1*N2,1)));
    rho_ic  = exp((x_ic-Vext)/kBT);
    figure('Color','white'); 
    subplot(1,2,1);
    HS.doPlots(rho_ic,'SC');
    subplot(1,2,2);
    opts.nContours = {0.1,0.2,0.3,0.4,0.5,0.6,0.7};
    HS.doPlots(rho_ic,'contour',opts.nContours);
    t_eqSol = toc;
    fprintf([num2str(t_eqSol),'s']);

    %************************************************
    %****************  Postprocess  ****************
    %************************************************        
    %SaveToFile(optsNum.DDFTCode,v2struct(data,optsPhys,optsNum),getResultsPath());    
    display(['Preprocessor, Computation time (sec): ', num2str(t_preprocess)]);
    display(['Equilibrium, Computation time (sec): ', num2str(t_eqSol)]);        
        
        
    %***************************************************************
    %   Physical Auxiliary functions:
    %***************************************************************         
    function x_ic = ComputeEquilibriumCondition(params,hs)
        x_ic    = fsolve(@GetExcessChemPotential,x_ig(mark));
        %kBT*log(rhoGas_sat*ones(N1*N2,1)));
    end
    

    function mu_s = GetExcessChemPotential(xm)
        
        x(mark)      = xm;
        x(~mark)     = x_ig(~mark);
        x            = x';
        
        rho_s  = exp((x-Vext)/kBT);        
        
        mu_s         = x + Conv*rho_s - mu + HS_f(rho_s,kBT) + Vadd;                
        mu_s         = mu_s(mark);
    end   

    function y = f_1DLG(x)
        y = [GetExcessChemPotential1D(x);x(ceil(end/2))-x_ig(ceil(end/2))];
    end

    function [om_ex,f_loc] = GetExcessGrandPotential(rho)
        %Compute excess grand potential
        f_id       = kBT*rho.*(log(rho)-1);
        [h1s,f_hs]   = HS_f(rho,kBT);
        f_attr     = 0.5*rho.*(Conv1D*rho);
        f_Vmu      = (Vext1D + Vadd1D - mu).*rho;        
        
        f_loc  = f_id + f_hs + f_attr + f_Vmu;                 
        om_ex  = Int_1D*(f_loc-f_loc(end));        
        
        %Get Bulk value and compare:
        floc_Bulk = GetfBulk(rho(end));
        fprintf(['Bulk Free Energy at inf : ',num2str(floc_Bulk),' with error: ',num2str(abs(floc_Bulk-f_loc(end))),'\n']);        
        floc_Bulk = GetfBulk(rho(1));
        fprintf(['Bulk Free Energy at -inf: ',num2str(floc_Bulk),' with error: ',num2str(abs(floc_Bulk-f_loc(1))),'\n']);         
        f_loc  = f_loc-f_loc(end);
    end

    function om_Bulk = GetfBulk(rho_Bulk)        
        [~,f_hs_Bulk]     = HS_f(rho_Bulk,kBT);
        Phi_r             = str2func(optsPhys.V2.V2DV2);        
        [h1s,h2s,alpha]       = Phi_r(0);    
        f_attr_Bulk       = alpha*rho_Bulk^2;
        f_id_Bulk         = kBT*rho_Bulk.*(log(rho_Bulk)-1);
        f_Vmu_Bulk        = -mu*rho_Bulk;
        
        om_Bulk           = f_id_Bulk + f_attr_Bulk + f_hs_Bulk + f_Vmu_Bulk;
    end

    function mu_s = GetExcessChemPotential1D(x)
        rho_s = exp((x-Vext1D)/kBT);
        mu_s  = x + Conv1D*rho_s - mu + HS_f(rho_s,kBT) + Vadd1D;
    end   

    function z = Phi(r)
         z = Phi_r(r,optsPhys.V2);
    end       

end