function PlotDensitySlices(this)       

    %% Initialization 
    global dirData
    
    n     = 5;
    y2Max = 15.5;
    y1Min = 0;
    y1Max = 15;
    
    rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;     
    
    y1P = y1Min + (y1Max-y1Min)*(0:1:n)/n;
    this.optsNum.PlotArea = struct('y1Min',y1Min,'y1Max',y1Max,...
                                   'y2Min',0.5,'y2Max',y2Max,...
                                   'N1',100,'N2',100);
    InitInterpolation(this);

    str = {'r','b','m','g','c'};
    
    %% Plotting
    f1 = figure('color','white','Position',[0 0 700 700]);    

    for i = 1:n
        % get adsorption

        ell      = this.HS.doIntFLine([y1P(i) y1P(i)],[0.5 y2Max],this.rho_eq-rhoGas_sat,'CHEB')/(rhoLiq_sat-rhoGas_sat);
        [rho,mu] = GetPointAdsorptionIsotherm(this,ell);        
        
        hold on;
        plot(this.AdsorptionIsotherm_Pts.y2_kv-0.5,rho,[str{i},'--'],'linewidth',1.5); hold on;
        this.HS.doPlotFLine([y1P(i) y1P(i)],[0.5 y2Max],this.rho_eq,'CART',str{i});        
    end
    
    box on;
    xlim([0 (y2Max-0.5)]);
    ylim([0 1.1]);%max(this.rho_eq)]);
    xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$n\sigma^3$','Interpreter','Latex','fontsize',25);
    set(gca,'fontsize',20);                        
    set(gca,'linewidth',1.5);          
    pbaspect([(y1Max-y1Min) (y2Max) 1]);

    f2 = figure('color','white','Position',[0 0 600 500]);    
    
    PlotEquilibriumResults(this,[y1Min y1Max],[0.5 y2Max],true,false); hold on;        
    for i = 1:n
        plot([y1P(i) y1P(i)],[0 (y2Max-0.5)],[str{i},':'],'linewidth',1.5);
    end

    inset2(f1,f2,0.4,[0.55,0.55]);
    set(gca,'fontsize',10);
    
    %inset2(f1,f2,0.35,[0.22,0.55]);
    close(f2);      
        
    %% Save plot
	print2eps([dirData filesep 'DensitySlices'],f1);
    saveas(f1,[dirData filesep 'DensitySlices.fig']);

end