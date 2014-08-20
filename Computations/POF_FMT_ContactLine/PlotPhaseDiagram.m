function PlotPhaseDiagram()  

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

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'maxComp_y2',10,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.49);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      
                  
    optsPhys.HSBulk = (['FexBulk_',optsNum.FexNum.Fex]);      

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
    %[rhoGas_satP,rhoLiq_satP,mu_satP,kBT_crit,rho_crit,mu_crit] = BulkPhaseDiagram(optsPhys);

    close all;
    f1 = figure('color','white','Position',[0 0 600 600]);
    
    plot(rhoGas_sat,kBT,'k','linewidth',1.5);hold on;
    plot(rhoLiq_sat,kBT,'k','linewidth',1.5);
    plot(rho_crit,kBT_crit,'ob','MarkerFace','b','MarkerSize',10);   
    if(isfield(optsPhys,'kBT'))
        plot(rhoGas_satP,optsPhys.kBT,'^r','MarkerFace','r','MarkerSize',10);
        plot(rhoLiq_satP,optsPhys.kBT,'^r','MarkerFace','r','MarkerSize',10);
    end
    xlabel('$n\sigma^3$','Interpreter','Latex','fontsize',fontS + 5); %\{\text{liq},\text{gas}\},\text{sat} 
    ylabel('$k_BT/\varepsilon$','Interpreter','Latex','fontsize',20);
    ylim([min(kBT),1.35])
    xlim([0 0.82]);
    
    set(gca,'ytick',0.5:0.1:1.4);
    set(gca,'fontsize',fontS);                        
    set(gca,'linewidth',1.5);  

        
    fid = fopen('D://Data/MichelsNGasNLiq.txt');
    x = textscan(fid,'%f %f %f','headerlines',3); %[T, rhoG, rhoL]
    fclose(fid); 
    plot(x{2},x{1},'kd','MarkerFaceColor','k','MarkerSize',10);
    plot(x{3},x{1},'kd','MarkerFaceColor','k','MarkerSize',10);
    hold on;	
    
    fid = fopen('D://Data/TrokhymchukNGas.txt');
    xGas = textscan(fid,'%f %f %f','headerlines',4); %[T, rhoG, rhoL]
    fclose(fid);
    
    fid = fopen('D://Data/TrokhymchukNLiq.txt');
    xLiq = textscan(fid,'%f %f %f','headerlines',3); %[T, rhoG, rhoL]
    fclose(fid);    
    plot(xGas{3},xGas{1},'ks','MarkerFaceColor','k','MarkerSize',10);
    plot(xLiq{3},xLiq{1},'ks','MarkerFaceColor','k','MarkerSize',10);
    hold on;
    
    print2eps([dirData filesep 'BulkPhaseDiagram'],f1);
    saveas(f1,[dirData filesep 'BulkPhaseDiagram.fig']);    
    
%            
%     %2 Trokhymchuk,  J.Chem.Phys. Vol. 111 No. 18. p.8510
%     fid = fopen('Data/TrokhymchukPressure.txt');
%     x = textscan(fid,'%f %f %f','headerlines',4); %[T, rhoG, rhoL]
%     fclose(fid);
%     subplot(2,2,4) 
%     plot(x{2},x{3},'rs');
%     hold on;
% 
%     
%     fid = fopen('Data/MichelsPressure.txt');
%     x = textscan(fid,'%f %f','headerlines',4); %[T, p]
     %fclose(fid);
%     subplot(2,2,4) 
%     plot(x{1},x{2},'o b');
%     hold on;
    

end