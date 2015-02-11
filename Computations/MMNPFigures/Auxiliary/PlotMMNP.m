function PlotMMNP
    close all;
        
    noRows  = 3;
    figMain = figure('color','white','Position',[0 0 1000 1300]);
    %s1 = subplot(4,2,2*(row-1)+1);
    
    SetRow('135deg',1,-10,'(a.I)','(a.II)',0.55,134.2);
    SetRow('120deg',2,-10,'(b.I)','(b.II)',0.7,119.9);
    SetRow('90deg',3,-10,'(c.I)','(c.II)',1.0,89.66);
    
    print2eps('D:\2DChebData\POF_FMT_ContactLine\DensitySlices_MMNPComparison\DensitySlices_Figure1',gcf);
    saveas(gcf,'D:\2DChebData\POF_FMT_ContactLine\DensitySlices_MMNPComparison\DensitySlices_Figure1.fig');   
    
    noRows  = 2;
    figMain = figure('color','white','Position',[0 0 1000 800]);
    SetRow('60deg',1,-10,'(d.I)','(d.II)',1.25,59.95);
    SetRow('40deg',2,-0,'(e.I)','(e.II)',1.375,41.02);             
    
    
    print2eps('D:\2DChebData\POF_FMT_ContactLine\DensitySlices_MMNPComparison\DensitySlices_Figure2',gcf);
    saveas(gcf,'D:\2DChebData\POF_FMT_ContactLine\DensitySlices_MMNPComparison\DensitySlices_Figure2.fig');   

    function SetRow(angle,row,xmin,str1,str2,alphaW,thetaY)
                
        f11 = openfig(['D:\2DChebData\POF_FMT_ContactLine\DensitySlices_MMNPComparison\',angle,'\DensitySlices_contour.fig']);
        ax11 = gca;                
        
        f12 = openfig(['D:\2DChebData\POF_FMT_ContactLine\DensitySlices_MMNPComparison\',angle,'\DensitySlices.fig']); %hgload
        ax12 = gca;        
                                
        
        f_cont = get(ax11,'children');
        
        set(0,'CurrentFigure',figMain)
        s1 = subplot(noRows,2,2*(row-1)+1); 
        copyobj(f_cont,s1);
        close(f11);
        
        xlim([xmin (xmin + 20)]);
        ylim([0 15]);

        xlabel('$x/\sigma$','Interpreter','Latex','fontsize',15);
        ylabel('$y/\sigma$','Interpreter','Latex','fontsize',15);
        box on;
        set(gca,'fontsize',10);
        pbaspect([20 15 1]);

        set(0,'CurrentFigure',figMain)
        s2 = subplot(noRows,2,2*(row-1)+2);        
        f_slices = get(ax12,'children');
        copyobj(f_slices,s2);        
        close(f12);
        xlabel('$y/\sigma$','Interpreter','Latex','fontsize',15);
        ylabel('$n\sigma^3$','Interpreter','Latex','fontsize',15);
        box on;
        set(gca,'fontsize',10);
        xlim([0 15]);
        ylim([0 1.4]);
        
        pbaspect([20 15 1]);
        
        
        ht1 = annotation('textbox',[0.07 (-0.05 + 0.9*(1-(row-1)/noRows)) 0.1 0.1]);
        set(ht1,'String',str1,'interpreter','none','FontSize',13);
        set(ht1','LineStyle','none');
        
        if(nargin >= 6)
            ht1 = annotation('textbox',[0.0 (-0.13 + 0.9*(1-(row-1)/noRows)) 0.1 0.1]);
            set(ht1,'String',['${\alpha_w \sigma^3}/{\varepsilon} = ',num2str(alphaW),'$'],'interpreter','Latex','FontSize',13);
            set(ht1','LineStyle','none');        
        end
        
        if(nargin >= 6)
            ht1 = annotation('textbox',[0.01 (-0.09 + 0.9*(1-(row-1)/noRows)) 0.1 0.1]);
            set(ht1,'String',['$\theta_Y = ',num2str(thetaY),'^\circ$'],'interpreter','Latex','FontSize',13);
            set(ht1','LineStyle','none');        
        end        
    
        ht1 = annotation('textbox',[0.5 (-0.05 + 0.9*(1-(row-1)/noRows)) 0.1 0.1]);
        set(ht1,'String',str2,'interpreter','none','FontSize',13);
        set(ht1','LineStyle','none'); 
	end
end