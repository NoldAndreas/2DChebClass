function FittingAdsorptionIsotherm(this,FT_Int,n)

    global dirData
    
    absFT = abs(this.AdsorptionIsotherm.FT);

    LogPlotInt = [2 FT_Int(2)];
    
    mark = ((absFT <= FT_Int(2)) & (absFT >= FT_Int(1))); %choose range
    
    ft = absFT(mark);
    mu = this.AdsorptionIsotherm.mu(mark)- this.optsPhys.mu_sat;    
            
    X = zeros(length(ft),n);
    
    for i = 1:n
        X(:,i) = ft.^(-(i+2));
    end
        
    a = (X'*X)\(X'*mu);
    this.AdsorptionIsotherm.Coeff = a;        
    disp(['coefficient = ',num2str(a)]);
    
    fft = LogPlotInt(1) + (0:0.01:1)'*(LogPlotInt(2)-LogPlotInt(1));
    fftX = zeros(length(fft),n);
    for i=1:n
        fftX(:,i) = fft.^(-(i+2));
    end
    
    mufft  = fftX*a;
    
    %Plot result
    f1 = figure('Color','white','Position',[0 0 800 800]);
    
    plot(this.AdsorptionIsotherm.mu- this.optsPhys.mu_sat,...
                    absFT,'k','linewidth',1.5); hold on;    
	plot(this.AdsorptionIsotherm.dmuCheck,...
                    absFT,'k:','linewidth',1.5); hold on;    
	plot([0 0],[0 max(absFT)],'k--','linewidth',1.5);
    set(gca,'fontsize',20);
    xlabel('$\Delta \mu/\varepsilon$','Interpreter','Latex','fontsize',25);
    ylabel('$\ell/\sigma$','Interpreter','Latex','fontsize',25);    
    ylim([0 max(absFT)]);    
    ylim([0 15]);
    xlim([min(this.AdsorptionIsotherm.dmuCheck)-0.02,max(this.AdsorptionIsotherm.dmuCheck)]);
    %xlim([-0.02 0.1]);
    
    f2 = figure('Color','white','Position',[0 0 800 800]); 
    set(gca,'YTickMode','manual')
    y = (abs(this.AdsorptionIsotherm.mu - this.optsPhys.mu_sat));
    
    markk = (absFT < LogPlotInt(2)) & (absFT > LogPlotInt(1)); 
    if(sum(markk)==0)
        return;
    end
    loglog((absFT(markk)),y(markk),'ko');  hold on;
        
    loglog(absFT(markk),...
                abs(this.AdsorptionIsotherm.dmuCheck(markk)),'k:','linewidth',1.5); hold on;    	
    
    hold on;
    loglog((fft),(abs(mufft)),'k--','linewidth',1.5);
    
    xlim(LogPlotInt);
    ylim([min(y(markk)),max(y(markk))]);
    ylim([min(y(markk)),0.1]);
    ylabel('$|\Delta \mu|/\varepsilon$','Interpreter','Latex','fontsize',25);
    xlabel('$\ell/\sigma$','Interpreter','Latex','fontsize',25);
    
    set(gca,'XTick',[2 5 10 15]); %LogPlotInt(1):5:LogPlotInt(2)
    set(gca,'YTick',[0.001 0.01 0.1]);%LogPlotInt(1):5:LogPlotInt(2));
    set(gca, 'XTickLabel', num2str(get(gca, 'XTick').'));
    set(gca, 'YTickLabel', num2str(get(gca, 'YTick').'));
    set(gca,'fontsize',20);
    
    if(min(this.AdsorptionIsotherm.FT)<0)
        inset2(f1,f2,0.45,[0.35,0.5]);
    else
        inset2(f1,f2,0.45,[0.5,0.5]);
    end
    close(f2);            
        
    folder = [dirData filesep 'EquilibriumSolutions'];    
    if(~exist(folder,'dir'))
                disp('Folder not found. Creating new path..');            
                mkdir(folder);
	end
    
    print2eps([folder filesep this.FilenameEq 'AdsorptionIsotherm'],f1);
    saveas(f1,[folder filesep this.FilenameEq 'AdsorptionIsotherm.fig']);
  
    
end