function [rhoGas_satP,rhoLiq_satP,mu_satP,kBT_crit,rho_crit,mu_crit] = BulkPhaseDiagram(optsPhys,filename)        
    
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
    figure('color','white','Position', [0 0 1500 500]);	
    fontS = 20;    
    
    subplot(1,3,1); PlotRho_T();    
	str = ['$T_{crit}$ = ',num2str(kBT_crit,3),', $\rho_{crit}$ = ',num2str(rho_crit,3)];
	%h = text(rho_crit-0.05,kBT_crit+0.02,['$T_{crit}$ = ',num2str(kBT_crit,3),', $\rho_{crit}$ = ',num2str(rho_crit,3)]);  set(h,'Interpreter','Latex'); set(h,'fontsize',fontS);
	title(str,'Interpreter','Latex','fontsize',fontS);
    
    subplot(1,3,2); Plot_T_ChemPot(); 
	str = ['$T_{crit}$ = ',num2str(kBT_crit,3),', $\mu_{crit}$ = ',num2str(mu_crit,3)];    
	%h = text(kBT_crit-0.1,mu_crit,str,'HorizontalAlignment','right');  set(h,'Interpreter','Latex'); set(h,'fontsize',fontS);
    title(str,'Interpreter','Latex','fontsize',fontS);
    
    subplot(1,3,3); Plot_T_P();
	str = ['$T_{crit}$ = ',num2str(kBT_crit,3),', $\mu_{crit}$ = ',num2str(mu_crit,3)];
	%h = text(kBT_crit-0.1,mu_crit,str,'HorizontalAlignment','right');  set(h,'Interpreter','Latex'); set(h,'fontsize',fontS);
	title(str,'Interpreter','Latex','fontsize',fontS);

    
    %Save Data     
    
    if((nargin >= 2) && ~isempty(filename))        
        h = figure('color','white','Position', [0 0 800 800]);	
        PlotRho_T();    
        print2eps([filename '_Bulk_Rho_T'],gcf);        
        saveas(gcf,[filename '_Bulk_Rho_T.fig']);   
        close(h);
        
        h = figure('color','white','Position', [0 0 800 800]);	
        Plot_T_ChemPot();           
        print2eps([filename '_Bulk_T_Mu'],gcf);        
        saveas(gcf,[filename '_Bulk_T_Mu.fig']);
        close(h);
        
        h = figure('color','white','Position', [0 0 800 800]);	
        Plot_T_P();        
        print2eps([filename '_Bulk_T_P'],gcf);        
        saveas(gcf,[filename '_Bulk_T_P.fig']);
        close(h);
    end    
    
    function PlotRho_T()
        plot(rhoGas_sat,kBT,'k','linewidth',1.5);hold on;
        plot(rhoLiq_sat,kBT,'k','linewidth',1.5);
        plot(rho_crit,kBT_crit,'ok','MarkerFace','k','MarkerSize',10);   
        if(isfield(optsPhys,'kBT'))
            plot(rhoGas_satP,optsPhys.kBT,'ok','MarkerFace','k','MarkerSize',10);
            plot(rhoLiq_satP,optsPhys.kBT,'ok','MarkerFace','k','MarkerSize',10);
        end
        xlabel('$\rho_{sat,\{gas,liq\}}$','Interpreter','Latex','fontsize',fontS+5); %\{\text{liq},\text{gas}\},\text{sat} 
        ylabel('$k_BT$','Interpreter','Latex','fontsize',fontS+5);
        ylim([min(kBT),max(kBT)+0.06])
        pbaspect([1 1 1]);   
        set(gca,'ytick',0.5:0.1:1);
        set(gca,'fontsize',fontS);                        
        set(gca,'linewidth',1.5);  
    end

    function Plot_T_ChemPot()
        %Plot Chemical Potential    
        plot(kBT,mu_sat,'k','linewidth',1.5); hold on;
        if(isfield(optsPhys,'kBT'))
            plot(optsPhys.kBT,mu_satP,'ok','MarkerFace','k','MarkerSize',10);
        end
        plot(kBT_crit,mu_crit,'ok','MarkerFace','k','MarkerSize',10);
        xlabel('$k_BT$','Interpreter','Latex','fontsize',fontS);
        ylabel('$\mu_{sat}$','Interpreter','Latex','fontsize',fontS);
        pbaspect([1 1 1]);   
        %set(gca,'ytick',-2.8:0.1:-2.5);
        %ylim([-2.85 -2.45]);
        set(gca,'fontsize',fontS);                        
        set(gca,'linewidth',1.5);       
    end

    function Plot_T_P()
        plot(kBT,p,'k','linewidth',1.5); hold on;
        if(isfield(optsPhys,'kBT'))
            plot(optsPhys.kBT,pP,'ok','MarkerFace','k','MarkerSize',10);
        end
        plot(kBT_crit,p_crit,'ok','MarkerFace','k','MarkerSize',10);
        xlabel('$k_BT$','Interpreter','Latex','fontsize',fontS);
        ylabel('$p_{sat}$','Interpreter','Latex','fontsize',fontS);
        pbaspect([1 1 1]);   
        %set(gca,'ytick',-2.8:0.1:-2.5);
        %ylim([-2.85 -2.45]);
        set(gca,'fontsize',fontS);                        
        set(gca,'linewidth',1.5);            
    end
end