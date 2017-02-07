function PlotDisjoiningPressureAnalysis(this)
    global dirData

    y1            = this.y1;
    filmThickness_II = this.IsolineInterface - 0.5;%filmThickness;
    filmThickness = this.filmThickness;
    f             = this.grandPot;    
    ST_1D         = this.ST_1D;
    
    theta_CS       = this.optsNum.PhysArea.alpha_deg*pi/180;
	rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;
      
	A      = 4*pi^2*(rhoLiq_sat-rhoGas_sat)*(rhoLiq_sat-this.optsPhys.V1.epsilon_w);
    fTAna  = (3:0.1:30)';
    wAna   = -A./(12*pi*(fTAna).^2);       
    dPAna  = A./(6*pi*(fTAna).^3);
    
    BulkPhaseDiagram(this.optsPhys);%,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq]);    
 
    %********************************************************************
    %Adsorption Isotherm based on 1D computations
    mu_sat = this.optsPhys.mu_sat;
    ell    = this.AdsorptionIsotherm_FT; 
    mu     = this.AdsorptionIsotherm_Mu;
    %rho    = this.AdsorptionIsotherm_rho;
	OmEx   = this.AdsorptionIsotherm_OmEx;
	dmuCheck = this.AdsorptionIsotherm_dmuCheck;

    figure('Color','white','Position',[0 0 800 800],'name','Adsorption Diagram');     
    plot([0 0],[min(ell) max(ell)],'k--','linewidth',1.5); hold on;
    plot(mu-mu_sat,ell,'k','linewidth',1.5); hold on;        
    plot(dmuCheck,ell,'k--','linewidth',1.5); hold on;       
    xlabel('$\Delta \mu$','Interpreter','Latex','fontsize',25);
    ylabel('$\ell$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20); set(gca,'linewidth',1.5);
    xlim([(min(mu- mu_sat)-0.01) (max(mu - mu_sat) + 0.01)]);
    ylim([(min(ell)-0.1) (max(ell)+0.1)]);

    print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_AdsorptionIsotherm'],gcf);
    saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_AdsorptionIsotherm.fig']);

    figure('Color','white','Position',[0 0 800 800]);     
    plot([mu_sat mu_sat],[min(OmEx) max(OmEx)],'r--','linewidth',1.5); hold on;
    plot(ell,OmEx,'linewidth',1.5); hold on;        
    title('Adsorption Diagram');
    xlabel('$\ell$','Interpreter','Latex','fontsize',25);
    ylabel('$F_{ex}$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20);
    set(gca,'linewidth',1.5);
    xlim([(min(ell)-0.01) (max(ell) + 0.01)]);
    ylim([(min(OmEx)-0.1) (max(OmEx)+0.1)]); 

    print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_AdsorptionIsothermOmega'],gcf);
    saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_AdsorptionIsothermOmega.fig']);
    
    %********************************************************************
    %1D Disjoining Pressure	    
    
    figure('Color','white','Position',[0 0 800 800],'name','Potentials');        
	plot(filmThickness,this.disjoiningPressure,'o'); hold on;
    plot([0,max(filmThickness)],[0,0],'k--');

    %Plot analytical solution including hamaker const:
    %plot(fTAna,dPAna,'r--','linewidth',1.5);        
    
    %Plot Comparison with Isotherm
    planarDisjoiningPressure = -(this.AdsorptionIsotherm_Mu-this.optsPhys.mu_sat)*(rhoLiq_sat-rhoGas_sat);
    
    plot(this.AdsorptionIsotherm_FT,planarDisjoiningPressure,'g','linewidth',1.5);
  %  plot(this.AdsorptionIsotherm_FT,this.disjoiningPressureCheck,'g--','linewidth',1.5);
    plot(this.AdsorptionIsotherm_FT,-dmuCheck*(rhoLiq_sat-rhoGas_sat),'g-.','linewidth',1.5);

    %plot(filmThickness,-this.ST_1D.om_LiqGas*(D2*filmThickness)./((1+(D*filmThickness).^2).^1.5));
%    plot(filmThickness,-this.ST_1D.om_LiqGas*(D2*filmThickness)./((1+(D*filmThickness).^2).^1.5),'b--');
   % plot(filmThickness_II(mark_II),-this.ST_1D.om_LiqGas*(D2_II*filmThickness_II(mark_II))./((1+(D_II*filmThickness_II(mark_II)).^2).^1.5),'r')
    
    title('Disjoining Pressure');
  %  ylim([1.2*min(this.disjoiningPressure) 1.2*max(this.disjoiningPressure)]);
   % xlim([min(filmThickness) max(filmThickness)]);       
    %********************************************************************
    %Disjoining Pressure as a function of the coordinate parallel to the
    %wall
    D  = this.DiffY1;
    D2 = this.DiffYY1;
    
    figure('Color','white','Position',[0 0 800 800],'name','Potentials');        
    
    plot(this.y1,zeros(size(this.y1)),'k--'); hold on;
    plot(this.y1,this.disjoiningPressure,'k:','linewidth',1.5); 
    
    mark = ((this.AdsorptionIsotherm_FT-min(this.AdsorptionIsotherm_FT))<max(filmThickness-min(filmThickness)));
    IP = barychebevalMatrix(filmThickness-min(filmThickness),this.AdsorptionIsotherm_FT(mark)-min(this.AdsorptionIsotherm_FT));
    plot(IP*this.y1,planarDisjoiningPressure(mark),'k-.','linewidth',1.5);
    
    plot(this.y1,-this.ST_1D.om_LiqGas*(D2*filmThickness)./((1+(D*filmThickness).^2).^1.5),'k--','linewidth',1.5);
    
    xlim([min(this.y1) max(this.y1)]);
    set(gca,'fontsize',20);
    xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$\Pi \sigma^3/\varepsilon$','Interpreter','Latex','fontsize',25);
    set(gca,'YTick',-0.1:0.01:0)
    
    print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_DisjoiningPressureX'],gcf);
    saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_DisjoiningPressureX.fig']);
    
%     
%     %********************************************************************
%     figure('Color','white','Position',[0 0 1000 1000],'name','Potentials');        
%     subplot(2,2,1);
%     %y1 = f_Box1.Pts.y1;                       
%     plot(y1,f,'o');  hold on;        
%     plot(y1,this.grandPot2,'--r');  hold on;       
%     plot([min(y1),max(y1)],[ST_1D.om_wallGas,ST_1D.om_wallGas],'k--');
%     %om  = ST_1D.om_wallLiq+ ST_1D.om_LiqGas/cos(alphaM);
%     %plot([min(f_Box1.Pts.y1),max(f_Box1.Pts.y1)],[om,om],'k--');        
%     om  = ST_1D.om_wallLiq+ ST_1D.om_LiqGas/cos(theta_CS);
%     plot([min(y1),max(y1)],[om,om],'g--');        
%     om  = ST_1D.om_wallLiq+ ST_1D.om_LiqGas/cos(this.alpha_YCA);
%     plot([min(y1),max(y1)],[om,om],'m--');        
%     title('Grand Potential(y1)');
%     xlabel('$y_1$','Interpreter','Latex');
% 
% %     subplot(2,2,2);    
% %     plot(filmThickness,f,'o'); hold on;       
% %     plot([0,max(filmThickness)],[ST_1D.om_wallGas,ST_1D.om_wallGas],'k--');
% %     om  = ST_1D.om_wallLiq+ ST_1D.om_LiqGas/cos(theta_CS);
% %     plot([min(filmThickness),max(filmThickness)],[om,om],'g--');        
% %     om  = ST_1D.om_wallLiq+ ST_1D.om_LiqGas/cos(this.alpha_YCA);
% %     plot([min(filmThickness),max(filmThickness)],[om,om],'m--');                       
% % %    om  = ST_1D.om_wallLiq+ ST_1D.om_LiqGas/cos(alphaM);
% % %    plot([0,max(hF)],[om,om],'k--');
% % 
% %     title('Grand Potential(h)');
% %     xlabel('Film Thickness');
% 
%     subplot(2,2,2);
%     plot(y1,this.disjoiningPressure,'o'); hold on;
%     plot([0,max(y1)],[0,0],'k--');   
%     %plot(fTAna,dPAna,'r--','linewidth',1.5);        
%     title('Disjoining Potential(y1)');
%     %ylim([1.2*min(w) 1.2*max(w)]);
%   %  xlim([0 10]);
% 
% 
%     %***************************************************************
%     %Compute Disjoining Potential
%     subplot(2,2,3);
%     
%     %m      = length(this.y1);
%     %e1     = [1,0.5*ones(1,m-2)];
%     %e2     = [0.5*ones(1,m-2),1];
%     %D      = diag(e1,1)-diag(e2,-1);
%     %D(1,1) = -1;  
%     %D(m,m) = 1;
%     %D      = D*(m-1)/(max(y1)-min(y1));    
%     %hF_P     = D*filmThickness;
%     
% 
% 
%     %Print Slope Function
%     subplot(2,2,4);
%     plot(y1,atan(hF_P)*180/pi,'o'); hold on;
%     plot(y1,ones(size(y1))*theta_CS*180/pi,'g--');
%     plot(y1,ones(size(y1))*this.alpha_YCA*180/pi,'m--');
%     
%     %Analytical Prediction of slope:
%    % slope = 180/pi*(acos( -wAna/this.ST_1D.om_LiqGas + cos(this.alpha_YCA) ));
%    % plot(y1,slope,'k-.','linewidth',1.5);
%     
%     title('Slope');
%     xlabel('$y_1$','Interpreter','Latex');
%     ylabel('[deg]');
%     ylim([1.1*min(atan(hF_P)*180/pi) 1.1*max(atan(hF_P)*180/pi)]);
% 
%     if(this.saveFigs)
%         print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_Potentials'],gcf);
%         saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_Potentials.fig']);
%     end  
    
    
end
    %****************************************************
    %****************************************************
%         [h1h,h2h,Int2] = HS.ComputeIntegrationVector();
%         Int2_Y = kron(eye(N1),Int2);
%         
%         [h1h,h2h,Int21_AD] = HS.AD.Sub_Strip.ComputeIntegrationVector();
%         [h1h,h2h,Int22_AD] = HS.AD.Sub_HalfSpace.ComputeIntegrationVector();
%                 
%         [h1h,h2h,Int2_AD] = HS.AD.ComputeIntegrationVector();
%         Int2_AD_Y         = [kron(eye(HS.AD.Sub_Strip.N1),Int21_AD) , ...
%                              kron(eye(HS.AD.Sub_HalfSpace.N1),Int22_AD)];
%                 
%         f_hs_end = [kron(f_hs(HS.AD.Pts.y2_kv==inf),ones(HS.AD.Sub_Strip.N2,1));...
%                     kron(f_hs(HS.AD.Pts.y2_kv==inf),ones(HS.AD.Sub_HalfSpace.N2,1))];
%         Fex    = Int2_AD_Y*(f_hs-f_hs_end)...
%                + Int2_Y*(f_loc-kron(f_loc(HS.Pts.y2_kv==inf),ones(N2,1))) +R*f_hs(Pts.y2_kv==inf);  

%    markY1     = (f_Box1.Pts.y2_kv==max(f_Box1.Pts.y2_kv));
%   markY3     = (f_Box3.Pts.y2_kv==max(f_Box3.Pts.y2_kv));
%    f_maxY2    = IP_1(markY1,:)*f_loc + IP_3(markY3,:)*f_hs;

%    PrintErrorPos((floc_Bulk+f_hs_Bulk)-f_maxY2,'Error of bulk free energy');%['Bulk Free Energy in 1D Computation for rho=',num2str(rho(end),2)]);