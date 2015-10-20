function PlotDensitySlices(this)       
    % Initialization         
    n     = 5;
    
    y2Max = this.optsNum.PlotAreaCart.y2Max;
    y1Min = this.optsNum.PlotAreaCart.y1Min;
    y1Max = this.optsNum.PlotAreaCart.y1Max;
        
    y1P = y1Min + (y1Max-y1Min)*(0.5:1:(n-0.5))/n;
    col = distinguishable_colors_NoRedBlueGreen();
    %str = {'','','m','','c'};
        
    % Plotting
    f1 = figure('color','white','Position',[0 0 700 700]);    
    

    for i = 1:n
        this.IDC.plotLine([y1P(i) y1P(i)],[0.5 y2Max],this.GetRhoEq,struct('dist0',true,'plain',true,'CART',true,'color',col(i,:)));  hold on;
    end
        
        % get adsorption
	if(this.optsNum.PhysArea.alpha_deg ~= 90)
        for i = 1:n
            ell      = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1P(i))*this.hIII;        
            [rho,mu] = GetPointAdsorptionIsotherm(this,ell);        
            mark     = (this.AdsorptionIsotherm.Pts.y2 < (y2Max + 0.5));
            plot(this.AdsorptionIsotherm.Pts.y2(mark)-0.5,rho(mark),'--','color',col(i,:)); hold on;%%%%
        end
    end
    
    box on;
    xlim([0 (y2Max-0.5)]);
    ylim([0 1.4]);%max(this.rho_eq)]);
    xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$n\sigma^3$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20);                        
    set(gca,'linewidth',1.5);          
    %pbaspect([(y1Max-y1Min) (y2Max) 1]);
    
    %SaveFigure('DensitySlices');
    SaveCurrentFigure(this,'DensitySlices');
    
    f2 = figure('color','white','Position',[0 0 600 500]);    
    if(this.optsNum.PhysArea.alpha_deg ~= 90)
        PlotContourResults(this,{'hI','hII','hIII'}); hold on;        
    else
        PlotContourResults(this,{'hI','hII'}); hold on;        
    end
    for i = 1:n
        plot([y1P(i) y1P(i)],[0 (y2Max-0.5)],':','color',col(i,:),'linewidth',1.3);
    end

    %SaveFigure('DensitySlices_contour');    
    SaveCurrentFigure(this,'DensitySlices_contour');    
    
    inset2(f1,f2,0.43,[0.55,0.62]);
    set(gca,'fontsize',10);
    
    %inset2(f1,f2,0.35,[0.22,0.55]);
    close(f2);      
        
    % Save plot
    SaveCurrentFigure(this,'DensitySlices_Inset');	
    %SaveCurrentFigure('DensitySlices_Inset');	
    
    function str = distinguishable_colors_NoRedBlueGreen()
        str = [0 1 1;... %cyan
               0.9412 0.4706 0;...	%Orange
               0.251 0 0.502;...%	Purple
               0 0.502 0.502;...%	Turquoise               
               1 0 1;...%pink
               1 0.502 0.502;...%peach	
               0.502 0.251 0]; %	Brown
    end

end