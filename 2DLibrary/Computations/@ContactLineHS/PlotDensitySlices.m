function PlotDensitySlices(this)       
    % Initialization 
    global dirData
    
    n     = 5;
    
    y2Max = this.optsNum.PlotAreaCart.y2Max;
    y1Min = this.optsNum.PlotAreaCart.y1Min;
    y1Max = this.optsNum.PlotAreaCart.y1Max;
        
    y1P = y1Min + (y1Max-y1Min)*(0:1:n)/n;
    str = {'r','b','m','g','c'};
    
    % Plotting
    f1 = figure('color','white','Position',[0 0 700 700]);    

    for i = 1:n
        % get adsorption
        ell      = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1P(i))*this.hIII;
        %ell      = this.IDC.doIntFLine([y1P(i) y1P(i)],[0.5 y2Max],this.GetRhoEq-rhoGas_sat,'CHEB')/(rhoLiq_sat-rhoGas_sat);
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);        
        
        hold on;
        plot(this.AdsorptionIsotherm.Pts.y2_kv-0.5,rho,[str{i},'--'],'linewidth',1.5); hold on;
        this.IDC.doPlotFLine([y1P(i) y1P(i)],[0.5 y2Max],this.GetRhoEq,'CART',str{i});        
    end
    
    box on;
    xlim([0 (y2Max-0.5)]);
    ylim([0 1.1]);%max(this.rho_eq)]);
    xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$n\sigma^3$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20);                        
    set(gca,'linewidth',1.5);          
    %pbaspect([(y1Max-y1Min) (y2Max) 1]);

    f2 = figure('color','white','Position',[0 0 600 500]);    
	PlotContourResults(this,true); hold on;        
    for i = 1:n
        plot([y1P(i) y1P(i)],[0 (y2Max-0.5)],[str{i},':'],'linewidth',1.5);
    end

    inset2(f1,f2,0.43,[0.55,0.62]);
    set(gca,'fontsize',10);
    
    %inset2(f1,f2,0.35,[0.22,0.55]);
    close(f2);      
        
    % Save plot
	print2eps([dirData filesep 'DensitySlices'],f1);
    saveas(f1,[dirData filesep 'DensitySlices.fig']);

end