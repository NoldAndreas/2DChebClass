function PlotDensitySlices(this)       
    % Initialization 
    global dirData
    
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
        % get adsorption
        ell      = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1P(i))*this.hIII;
        %ell      = this.IDC.doIntFLine([y1P(i) y1P(i)],[0.5 y2Max],this.GetRhoEq-rhoGas_sat,'CHEB')/(rhoLiq_sat-rhoGas_sat);
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);        
        
        hold on;
        plot(this.AdsorptionIsotherm.Pts.y2-0.5,rho,'--','color',col(i,:),'linewidth',1.5); hold on;%%%%
        this.IDC.plotLine([y1P(i) y1P(i)],[0.5 y2Max],this.GetRhoEq,struct('dist0',true,'plain',true,'CART',true,'color',col(i,:)));        
    end
    
    box on;
    xlim([0 (y2Max-0.5)]);
    ylim([0 1.4]);%max(this.rho_eq)]);
    xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$n\sigma^3$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20);                        
    set(gca,'linewidth',1.5);          
    %pbaspect([(y1Max-y1Min) (y2Max) 1]);
    
    SaveFigure('DensitySlices');
    
    f2 = figure('color','white','Position',[0 0 600 500]);    
	PlotContourResults(this,true); hold on;        
    for i = 1:n
        plot([y1P(i) y1P(i)],[0 (y2Max-0.5)],':','color',col(i,:),'linewidth',1.5);
    end

    print2eps([dirData filesep 'DensitySlices_contour'],f2);
    saveas(f2,[dirData filesep 'DensitySlices_contour.fig']);
    
    inset2(f1,f2,0.43,[0.55,0.62]);
    set(gca,'fontsize',10);
    
    %inset2(f1,f2,0.35,[0.22,0.55]);
    close(f2);      
        
    % Save plot
	print2eps([dirData filesep 'DensitySlices_Inset'],f1);
    saveas(f1,[dirData filesep 'DensitySlices_Inset.fig']);
    
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