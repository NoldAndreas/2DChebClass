function PhaseDiagrams()  

    fontS = 15;

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    

    PhysArea = struct('N',[45,90],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',20);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);

    %Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
%                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
    optsNum = struct('PhysArea',PhysArea,...                     
                     'maxComp_y2',10,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.49);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'HSBulk','CarnahanStarling',...
                      'sigmaS',1);      
                     

    %******************************************
    %******************************************    
    [kBT_crit,rho_crit,mu_crit,p_crit] = GetCriticalPoint(optsPhys,[],false);
        
    intitialGuess = [0.01;0.6;-2];
    i = 1;
    optsPhysVar     = optsPhys;
    optsPhysVar.kBT = kBT_crit/2;
    
    while((intitialGuess(1) ~= intitialGuess(2)) && (optsPhysVar.kBT < kBT_crit))
        
        [rhoGas_sat(i),rhoLiq_sat(i),mu_sat(i),p(i)] = BulkSatValues(optsPhysVar,intitialGuess,false);
        kBT(i)       = optsPhysVar.kBT;
        initialGuess = [rhoGas_sat(i),rhoLiq_sat(i),mu_sat(i)];
        
        
        optsPhysVar.kBT = optsPhysVar.kBT + 0.005;
        i = i + 1;
    end
    rhoGas_sat = real(rhoGas_sat);
    rhoLiq_sat = real(rhoLiq_sat);
    mu_sat     = real(mu_sat);
    p          = real(p);
    
    rhoGas_sat(i) = rho_crit;
    rhoLiq_sat(i) = rho_crit;
    mu_sat(i)     = mu_crit;
    kBT(i)        = kBT_crit;
    p(i)          = p_crit;
    
    if(isfield(optsPhys,'kBT'))
        [rhoGas_satP,rhoLiq_satP,mu_satP,pP] = BulkSatValues(optsPhys,intitialGuess,false);
    end
    %******************************************
    %******************************************
    optsPhys.criticalValues = false;
    [rhoGas_satP,rhoLiq_satP,mu_satP,kBT_crit,rho_crit,mu_crit] = BulkPhaseDiagram(optsPhys);
    
    subplot(1,3,1);
	xlabel('$v_0 n$','Interpreter','Latex','fontsize',fontS + 5); %\{\text{liq},\text{gas}\},\text{sat} 
    ylabel('$k_BT/\varepsilon$','Interpreter','Latex','fontsize',20);
    
    subplot(1,3,2);
	xlabel('$v_0 n$','Interpreter','Latex','fontsize',fontS + 5); %\{\text{liq},\text{gas}\},\text{sat} 
    ylabel('$\mu_{sat}/\varepsilon$','Interpreter','Latex','fontsize',20);

    subplot(1,3,3);
	xlabel('$v_0 n$','Interpreter','Latex','fontsize',fontS + 5); %\{\text{liq},\text{gas}\},\text{sat} 
    ylabel('$v_0 p/\varepsilon$','Interpreter','Latex','fontsize',20);
    
	print2eps(['CarnahanStarling_PhaseDiagram'],gcf);        
	saveas(gcf,['CarnahanStarling_PhaseDiagram.fig']);    

end