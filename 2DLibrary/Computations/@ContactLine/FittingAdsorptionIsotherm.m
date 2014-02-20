function FittingAdsorptionIsotherm(this,FT_Int,n)

    global dirData

    LogPlotInt = [5 FT_Int(2)];
    
    mark = ((this.AdsorptionIsotherm_FT <= FT_Int(2)) & (this.AdsorptionIsotherm_FT >= FT_Int(1))); %choose range
    
    ft = this.AdsorptionIsotherm_FT(mark);
    mu = this.AdsorptionIsotherm_Mu(mark)- this.optsPhys.mu_sat;    
            
    X = zeros(length(ft),n);
    
    for i = 1:n
        X(:,i) = ft.^(-(i+2));
    end
        
    a = (X'*X)\(X'*mu);
    this.AdsorptionIsotherm_Coeff = a;
    a
    
    
    fft = LogPlotInt(1) + (0:0.01:1)'*(LogPlotInt(2)-LogPlotInt(1));
    fftX = zeros(length(fft),n);
    for i=1:n
        fftX(:,i) = fft.^(-(i+2));
    end
    
    mufft  = fftX*a;
    
    %Plot result
    f1 = figure('Color','white','Position',[0 0 800 800]);
    
    plot(this.AdsorptionIsotherm_Mu- this.optsPhys.mu_sat,...
                    this.AdsorptionIsotherm_FT,'k','linewidth',1.5); hold on;    
	plot(this.AdsorptionIsotherm_dmuCheck,...
                    this.AdsorptionIsotherm_FT,'k:','linewidth',1.5); hold on;    
	plot([0 0],[0 max(this.AdsorptionIsotherm_FT)],'k--','linewidth',1.5);
    set(gca,'fontsize',20);
    xlabel('$\Delta \mu/\varepsilon$','Interpreter','Latex','fontsize',25);
    ylabel('$\ell/\sigma$','Interpreter','Latex','fontsize',25);    
    ylim([0 max(this.AdsorptionIsotherm_FT)]);
    ylim([0 15]);
    xlim([-0.02 0.1]);
    
    f2 = figure('Color','white','Position',[0 0 800 800]); 
    set(gca,'YTickMode','manual')
    y = (abs(this.AdsorptionIsotherm_Mu - this.optsPhys.mu_sat));
    
    markk = (this.AdsorptionIsotherm_FT < LogPlotInt(2)) & (this.AdsorptionIsotherm_FT > LogPlotInt(1)); 
    loglog((this.AdsorptionIsotherm_FT(markk)),y(markk),'ko');  hold on;
        
    loglog(this.AdsorptionIsotherm_FT(markk),...
                abs(this.AdsorptionIsotherm_dmuCheck(markk)),'k:','linewidth',1.5); hold on;    	
    
    hold on;
    loglog((fft),(abs(mufft)),'k--','linewidth',1.5);
    
    xlim(LogPlotInt);
    ylim([min(y(markk)),max(y(markk))]);
    ylabel('$|\Delta \mu|/\varepsilon$','Interpreter','Latex','fontsize',25);
    xlabel('$\ell/\sigma$','Interpreter','Latex','fontsize',25);
    
    set(gca,'XTick',LogPlotInt(1):5:LogPlotInt(2));
    set(gca,'YTick',0.0005:0.0005:0.002);%LogPlotInt(1):5:LogPlotInt(2));
    set(gca, 'XTickLabel', num2str(get(gca, 'XTick').'));
    set(gca, 'YTickLabel', num2str(get(gca, 'YTick').'));
    set(gca,'fontsize',20);
    
    inset2(f1,f2,0.45,[0.5,0.5]);
    close(f2);            
    
    print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_AdsorptionIsotherm'],f1);
    saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_AdsorptionIsotherm.fig']);
  
    
end