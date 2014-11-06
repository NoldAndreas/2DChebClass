function f2 = PlotDisjoiningPressures(this)

	y1   = this.y1_SpectralLine.Pts.y;
    fT   = abs(this.AdsorptionIsotherm.FT);
    D    = this.y1_SpectralLine.Diff.Dy;
    D2   = this.y1_SpectralLine.Diff.DDy;
    
    Pi_I = GetDisjoiningPressure_I(this);
        
    f2 = figure('Color','white','Position',[0 0 1000 750]);            
    
    plot([-10 35],[0 0],'k'); hold on;
    plot(y1,this.disjoiningPressure_II,'k--','linewidth',1.5);     
    plot(y1,this.disjoiningPressure_IV,'r--','linewidth',1.5);     
    
    [DeltaY1_II,DeltaY1_III] = this.ComputeDeltaFit();
    plot(this.y1_I+DeltaY1_II,disjoiningPressure1D(this.hI),'k-.','linewidth',1.5);
    plot(this.y1_I+DeltaY1_III,disjoiningPressure1D(this.hI),'k-.','linewidth',1.5);
    
   % mark = ((this.AdsorptionIsotherm_FT-min(this.AdsorptionIsotherm_FT))<max(filmThickness-min(filmThickness)));
   % IP = barychebevalMatrix(filmThickness-min(filmThickness),this.AdsorptionIsotherm_FT(mark)-min(this.AdsorptionIsotherm_FT));
   % plot(IP*this.y1,planarDisjoiningPressure(mark),'k-.','linewidth',1.5);
   
    markFT = (this.hIII < 25);
    hhFT   = -this.ST_1D.om_LiqGas*(D2*this.hIII)./((1+(D*this.hIII).^2).^1.5);
   
    plot(y1(markFT),hhFT(markFT),'k','linewidth',1.5);
    
    xlim([min(y1) max(y1)]);
    %ylim([1.1*min(dp1D) 3*max(0.001,max(dp1D))]);%[-0.035 0.005])
    set(gca,'fontsize',20);
    xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$\Pi \sigma^3/\varepsilon$','Interpreter','Latex','fontsize',25);
    %set(gca,'YTick',-0.1:0.01:0);    
    
    global dirData
    print2eps([dirData filesep 'Equilibrium' filesep this.FilenameEq '_DisjoiningPressures'],gcf);
	saveas(gcf,[dirData filesep 'Equilibrium' filesep this.FilenameEq '_DisjoiningPressures.fig']);

    function dP1D = disjoiningPressure1D(ell) 
        
        ellAD = min(fT)+ell;
        IP    = barychebevalMatrix(fT,ellAD);
        dP1D  = IP*Pi_I;
        dP1D(ellAD<min(fT)) = 0;
        dP1D(ellAD>max(fT)) = 0;
    end
    
end