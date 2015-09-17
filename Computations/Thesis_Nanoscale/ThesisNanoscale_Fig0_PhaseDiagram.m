function ThesisNanoscale_Fig0_PhaseDiagram()  

    fontS = 15;
    
    AddPaths('ThesisNanoscale');
    
	config = ThesisNanoscale_GetStandardConfig([],[]);
    config.optsPhys.V2 = rmfield(config.optsPhys.V2,'r_cutoff');
 
    optsPhys        = config.optsPhys;              
    optsPhys.HSBulk = (['FexBulk_',config.optsNum.FexNum.Fex]);      

    %******************************************
    %******************************************    
    [kBT_crit,rho_crit,mu_crit,p_crit] = GetCriticalPoint(optsPhys,[]);
        
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
        [rhoGas_satP,rhoLiq_satP,mu_satP,pP] = BulkSatValues(optsPhys,intitialGuess);
    end
    %******************************************
    %******************************************
    [rhoGas_satP,rhoLiq_satP,mu_satP,kBT_crit,rho_crit,mu_crit] = BulkPhaseDiagram(optsPhys,{'PlotExperimentalData'});
        
    close all;
    f1 = figure('color','white','Position',[0 0 3500 300]);
    
    plot(rhoGas_sat,kBT,'k','linewidth',1.5);hold on;
    plot(rhoLiq_sat,kBT,'k','linewidth',1.5);
    plot(rho_crit,kBT_crit,'ob','MarkerFace','b','MarkerSize',10);   
    if(isfield(optsPhys,'kBT'))
        plot(rhoGas_satP,optsPhys.kBT,'^r','MarkerFace','r');%,'MarkerSize',10);
        plot(rhoLiq_satP,optsPhys.kBT,'^r','MarkerFace','r');%,'MarkerSize',10);
    end
    xlabel('$\nDensity$','Interpreter','Latex'); %\{\text{liq},\text{gas}\},\text{sat} 
    ylabel('$T$','Interpreter','Latex');
    ylim([min(kBT),1.35])
    xlim([0 0.82]);
    
    set(gca,'ytick',0.5:0.1:1.4);
    set(gca,'fontsize',fontS);                        
    set(gca,'linewidth',1.5);  

    PlotData('D://Data/MichelsNGasNLiq.txt','d');
    PlotData('D://Data/Trokhymchuk_cutoff_2_5.txt','s');
      
       
    SaveFigure('BulkPhaseDiagram');    

    
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