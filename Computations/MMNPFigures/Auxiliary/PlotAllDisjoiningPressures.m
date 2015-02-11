function PlotAllDisjoiningPressures()

        figMain = figure('color','white','Position',[0 0 1000 400]);
        row    = 1;
        noRows = 1;
        xmin = -10;

        f11 = openfig('D:\2DChebData\POF_FMT_ContactLine\DijoiningPressures_MMNPComparison\DisjoiningPressureMMNP_drying.fig');
        ax11 = gca;                
        
        f12 = openfig('D:\2DChebData\POF_FMT_ContactLine\DijoiningPressures_MMNPComparison\DisjoiningPressureMMNP_wetting.fig');
        ax12 = gca;        
                                
        
        f_cont = get(ax11,'children');
        
        set(0,'CurrentFigure',figMain)
        s1 = subplot(1,2,2*(row-1)+1); 
        copyobj(f_cont,s1);
        close(f11);
        
        xlim([xmin (xmin + 20)]);     
        ylim([-0.17 0.005])

        xlabel('$x/\sigma$','Interpreter','Latex','fontsize',20);
        ylabel('$\Pi\sigma^3/\varepsilon$','Interpreter','Latex','fontsize',20);
        box on;
        set(gca,'fontsize',15);

        set(0,'CurrentFigure',figMain)
        s2 = subplot(noRows,2,2*(row-1)+2);        
        f_slices = get(ax12,'children');
        copyobj(f_slices,s2);        
        close(f12);
        xlabel('$x/\sigma$','Interpreter','Latex','fontsize',20);
        ylabel('$\Pi\sigma^3/\varepsilon$','Interpreter','Latex','fontsize',20);
        box on;
        set(gca,'fontsize',15);
        xlim([-8 12]);        
        ylim([-0.17 0.005])                
        
        
        ht1 = annotation('textbox',[0.07 0.85 0.1 0.1]);
        set(ht1,'String','(a)','interpreter','none','FontSize',13);
        set(ht1','LineStyle','none');
    
        ht1 = annotation('textbox',[0.5 0.85 0.1 0.1]);
        set(ht1,'String','(b)','interpreter','none','FontSize',13);
        set(ht1','LineStyle','none'); 
        
        print2eps('D:\2DChebData\POF_FMT_ContactLine\DijoiningPressures_MMNPComparison\DisjoiningPressures',gcf);
        saveas(gcf,'D:\2DChebData\POF_FMT_ContactLine\DijoiningPressures_MMNPComparison\DisjoiningPressures.fig');   

        


end