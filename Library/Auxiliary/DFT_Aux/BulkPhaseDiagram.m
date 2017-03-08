function [rhoGas_satP,rhoLiq_satP,mu_satP,kBT_crit,rho_crit,mu_crit] = BulkPhaseDiagram(optsPhys,opts)
    if(nargin < 2)
        opts = {};
    end

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
    
    %*******************************************************************
    %Plot Density Profiles
    figure('color','white','Position', [0 0 250 200]);	                
    PlotRho_T();    
    if(IsOption(opts,'PlotExperimentalData'))
        PlotData('D://Data/MichelsNGasNLiq.txt','d'); hold on;
        PlotData('D://Data/Trokhymchuk_cutoff_5.txt','s');
        ylim([0.5 1.4]);
    end   
    SaveFigure('BulkPhaseDiagram_Bulk_Rho_T');
   
    figure('color','white','Position', [0 0 250 200]);	                
	Plot_T_ChemPot();           
    SaveFigure('BulkPhaseDiagram_Bulk_T_mu');        
        
    figure('color','white','Position', [0 0 250 200]);
	Plot_T_P();        
    SaveFigure('BulkPhaseDiagram_Bulk_T_p');
        
    function PlotRho_T()
        plot(rhoGas_sat,kBT,'k','linewidth',1.5);hold on;
        plot(rhoLiq_sat,kBT,'k','linewidth',1.5);
        plot(rho_crit,kBT_crit,'or','MarkerFace','r');   
        if(isfield(optsPhys,'kBT'))
            plot(rhoGas_satP,optsPhys.kBT,'ob','MarkerFace','b');
            plot(rhoLiq_satP,optsPhys.kBT,'ob','MarkerFace','b');
        end
        xlabel('$\nDensityV/\nDensityL$','Interpreter','Latex'); %\{\text{liq},\text{gas}\},\text{sat} 
        ylabel('$T$','Interpreter','Latex');
       % ylim([min(kBT),max(kBT)+0.06])
        pbaspect([1 1 1]);   
        %set(gca,'ytick',0.5:0.1:1);        
    end
    function Plot_T_ChemPot()
        %Plot Chemical Potential    
        plot(kBT,mu_sat,'k','linewidth',1.5); hold on;
        if(isfield(optsPhys,'kBT'))
            plot(optsPhys.kBT,mu_satP,'ob','MarkerFace','b');
        end
        plot(kBT_crit,mu_crit,'or','MarkerFace','r');
        xlabel('$T$','Interpreter','Latex');
        ylabel('$\chemPot_{\text{sat}}$','Interpreter','Latex');
        pbaspect([1 1 1]);   
        xlim([min(kBT),max(kBT)+0.06])
        %set(gca,'ytick',-2.8:0.1:-2.5);
        %ylim([-2.85 -2.45]);
    end
    function Plot_T_P()
        plot(kBT,p,'k','linewidth',1.5); hold on;
        if(isfield(optsPhys,'kBT'))
            plot(optsPhys.kBT,pP,'ob','MarkerFace','b');
        end
        plot(kBT_crit,p_crit,'or','MarkerFace','r');
        xlabel('$T$','Interpreter','Latex');
        ylabel('$p_{\text{sat}}$','Interpreter','Latex');
        pbaspect([1 1 1]);   
        xlim([min(kBT),max(kBT)+0.06])
        %set(gca,'ytick',-2.8:0.1:-2.5);
        %ylim([-2.85 -2.45]);        
    end

    function PlotData(filename,sym)
        fid = fopen(filename);
        x = textscan(fid,'%f %f %f','headerlines',3); %[T, rhoG, rhoL]
        fclose(fid); 
        plot(x{2},x{1},['k',sym],'MarkerFaceColor','k'); hold on;	
        plot(x{3},x{1},['k',sym],'MarkerFaceColor','k'); hold on;	       
    end
    

end