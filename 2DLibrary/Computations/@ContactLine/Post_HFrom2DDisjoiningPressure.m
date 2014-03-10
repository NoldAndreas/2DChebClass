function [f1,f2] = Post_HFrom2DDisjoiningPressure(this,f1)

    global dirData
    %*************************************
    %Initialization
    y1   = this.y1;
    DP   = this.disjoiningPressure;
    ST   = this.ST_1D.om_LiqGas;
    mu_sat         = this.optsPhys.mu_sat;    
    rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;
    planarDisjoiningPressure = -(this.AdsorptionIsotherm_Mu-mu_sat)*(rhoLiq_sat-rhoGas_sat);
    %*************************************
    
    %solves Eq. (17) from Henderson (2010)    
    %sets ell_2DDisjoiningPressure
    %uses disjoiningPressure,y1,Int_y1
    
    Diff = barychebdiff(y1,2);    
    
    %**Alternative Computation
	IntMY1 = zeros(length(y1));
	vh     = zeros(1,length(y1));
    for i = 2:length(y1)
        %Matrix for integral int(v(y),y=y1..Y)
        h           = y1(i) - y1(i-1);
        vh([i-1,i]) = vh([i-1,i]) + h/2;
        IntMY1(i,:)   = vh;
    end
    hP = -tan(asin(IntMY1*DP/ST));
    h_2D = fsolve(@f3,zeros(size(y1)));      
    this.ell_2DDisjoiningPressure   = h_2D;
    
    
    %*************************************
    %h-Profile from 1DDisjoiningPressure
	
    %****************************************
    %a) Compute Interface Potential from Disjoining Pressure
    
    %Film thickness   
    fT = this.AdsorptionIsotherm_FT;
	IntM = zeros(length(fT));
	vh   = zeros(1,length(fT));
    for i = 2:length(fT)        
        %Matrix for integral int(v(y),y=y1..Y)
        h           = fT(i) - fT(i-1);
        vh([i-1,i]) = vh([i-1,i]) + h/2;
        IntM(i,:)   = vh;
    end
    
    %Interface Potential
    IntPot = -IntM*planarDisjoiningPressure; 
    
   % plot(fT,planarDisjoiningPressure); hold on;
%    plot(fT,IntPot,'m');
    checkSumRule = this.ST_1D.om_LiqGas*(1-cos(this.alpha_YCA));
 %   plot([min(fT) max(fT)],[checkSumRule checkSumRule],'k--');
    
    x         = ClenCurtFlip(150);
    CompDiff  = barychebdiff(x,2);    
    
    if(this.optsNum.PhysArea.alpha_deg == 60)
        hmax =  this.filmThickness(75)-min(this.filmThickness);    
        L = this.y1(75)-min(this.y1);
        z = (x+1)/2*L -L+y1(75);
    else
        hmax =  this.filmThickness(end)-min(this.filmThickness);    
        L = max(this.y1)-min(this.y1);
        z = (x+1)/2*L -L+max(y1);
    end
    DPhys.Dx  = CompDiff.Dx*(2/L);    
    ST_LG = this.ST_1D.om_LiqGas;
    h  = fsolve(@f6,hmax*ones(size(z)));%exp(-z));%ell_h);
    
    %this.ell_1DDisjoiningPressure   = h; 
    
    
    hmax = h_2D(end);
    L2 = 10;
    z2 = (x+1)/2*L2 -L2 +max(this.y1);
    DPhys.Dx  = CompDiff.Dx*(2/L2);        
    h2   = fsolve(@f6,hmax*ones(size(z)));%exp(-z));%ell_h);
    
    hmax = h2(1);
    L3 = 12;
    z3 = (x+1)/2*L3 - L3 +min(z2);
    DPhys.Dx  = CompDiff.Dx*(2/L3);
    h3   = fsolve(@f6,hmax*ones(size(z)));%exp(-z));%ell_h);
    
    hmax = h3(1);
    L4 = 12;
    z4 = (x+1)/2*L4 - L4 +min(z3);
    DPhys.Dx  = CompDiff.Dx*(2/L4);
    h4   = fsolve(@f6,hmax*ones(size(z)));%exp(-z));%ell_h);
    
    
    %**********************
    %PLOT * Start
	%**********************
    %this.optsNum.PlotArea.y1Min = min(this.y1);
	%this.optsNum.PlotArea.y1Max = max(this.y1);
    this.optsNum.PlotArea.y1Min = -5;%min(this.y1);
	this.optsNum.PlotArea.y1Max = 30;%max(this.y1);
    this.optsNum.PlotArea.y2Min = 0.5;
    this.optsNum.PlotArea.y2Max = 25.5;%max(this.filmThickness)+2;
        
	this.optsNum.PlotArea.N1    = 130;
	this.optsNum.PlotArea.N2    = 130;
    
	InitInterpolation(this,true);
    
    rhoLiq_sat    = this.optsPhys.rhoLiq_sat;
	rhoGas_sat    = this.optsPhys.rhoGas_sat;
    h0 = min(this.filmThickness);
    
    ratio_y_x =  (this.optsNum.PlotArea.y2Max- this.optsNum.PlotArea.y2Min)/...
                    (this.optsNum.PlotArea.y1Max- this.optsNum.PlotArea.y1Min);
                
    LP_1 = this.optsNum.PlotArea.y1Max - this.optsNum.PlotArea.y1Min;
    LP_2 = this.optsNum.PlotArea.y2Max - this.optsNum.PlotArea.y2Min;
    LP   = 20;
    
	if(nargin < 2)
       % f1 = figure('Color','white','Position',[0 0 LP*(8+LP_1) LP*(8+LP_2)]);
        f1 = figure('Color','white','Position',[0 0 800 600]);
        %if(ratio_y_x < 1)
%            f1 = figure('Color','white','Position',[0 0 800/ratio_y_x 800]);
%        else
%            f1 = figure('Color','white','Position',[0 0 800 800*ratio_y_x]);
%        end
    end
    
    %subplot(2,1,1);
    plot(this.y1,h0 + h_2D,'k:','linewidth',1.5); hold on;
    
    plot(z,h0 + h,'k-.','linewidth',1.5);
    plot(z2,h0 + h2,'k-.','linewidth',1.5);
    plot(z3,h0 + h3,'k-.','linewidth',1.5);
    plot(z4,h0 + h4,'k-.','linewidth',1.5);
    %plot(z,h0 + h,'om');        
    plot(this.y1,this.filmThickness,'k','linewidth',1.5);            
    
    optDetails.clabel = false;  
	optDetails.linecolor = 'k';
    drho = rhoLiq_sat - rhoGas_sat;
    
    optDetails.y2CartShift = -0.5;
    optDetails.linewidth = 1;
	optDetails.nContours = rhoGas_sat + [0.05,0.5,0.95]*drho;
	%this.HS.doPlots(this.rho_eq,'contour',optDetails);  hold on;  
    PlotEquilibriumResults(this,[],[],true,false);
        
	xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
	ylabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    xlim([this.optsNum.PlotArea.y1Min this.optsNum.PlotArea.y1Max]);   
    %pbaspect([2 1 1]);
    
	%print2eps([dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_interfaceComparison'],gcf);
	%saveas(gcf,[dirData filesep 'EquilibriumSolutions' filesep this.FilenameEq '_interfaceComparison.fig']);
	%**********************
    %PLOT * End
	%**********************
    %**********************
    %PLOT DisjoiningPressure * Start
	%**********************
     
  %  ylim([1.2*min(this.disjoiningPressure) 1.2*max(this.disjoiningPressure)]);
   % xlim([min(filmThickness) max(filmThickness)]);       
    %********************************************************************
    %Disjoining Pressure as a function of the coordinate parallel to the
    %wall
    D = this.DiffY1;
    D2 = this.DiffYY1;
     
    %subplot(2,1,2);
    f2 = figure('Color','white','Position',[0 0 1000 750]);        
    
    dp1D = disjoiningPressure1D(h);
    plot(this.y1,zeros(size(this.y1)),'k--'); hold on;
    plot(this.y1,this.disjoiningPressure,'k:','linewidth',1.5); 
    plot(z,dp1D,'k-.','linewidth',1.5);
    plot(z2,disjoiningPressure1D(h2),'k-.','linewidth',1.5);
    plot(z3,disjoiningPressure1D(h3),'k-.','linewidth',1.5);
    plot(z4,disjoiningPressure1D(h4),'k-.','linewidth',1.5);
    
   % mark = ((this.AdsorptionIsotherm_FT-min(this.AdsorptionIsotherm_FT))<max(filmThickness-min(filmThickness)));
   % IP = barychebevalMatrix(filmThickness-min(filmThickness),this.AdsorptionIsotherm_FT(mark)-min(this.AdsorptionIsotherm_FT));
   % plot(IP*this.y1,planarDisjoiningPressure(mark),'k-.','linewidth',1.5);
   
    markFT = (this.filmThickness < 25);
    hhFT   = -this.ST_1D.om_LiqGas*(D2*this.filmThickness)./((1+(D*this.filmThickness).^2).^1.5);
   
    plot(this.y1(markFT),hhFT(markFT),'k','linewidth',1.5);
    
    xlim([min(this.y1) max(this.y1)]);
    ylim([1.1*min(dp1D) 3*max(0.001,max(dp1D))]);%[-0.035 0.005])
    set(gca,'fontsize',20);
    xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$\Pi \sigma^3/\varepsilon$','Interpreter','Latex','fontsize',25);
    %set(gca,'YTick',-0.1:0.01:0);
    %pbaspect([2 1 1]);      
    
    %**********************
    %PLOT DisjoiningPressure * End
	%**********************
    
    function y = f3(h)
        y    = Diff.Dx*h-hP;
        y(1) = h(1);
    end

    function y = f6(h)
        dh = ((1-interfacePotential(h)/ST_LG).^(-2)-1).^0.5;
        y  = DPhys.Dx*h - dh;
        y(end) = h(end) - hmax;
    end


    function dP1D = disjoiningPressure1D(ell) 
        
        ellAD = min(this.AdsorptionIsotherm_FT)+ell;
        IP    = barychebevalMatrix(this.AdsorptionIsotherm_FT,ellAD);
        dP1D  = IP*planarDisjoiningPressure;
        dP1D(ellAD<min(this.AdsorptionIsotherm_FT)) = 0;
        dP1D(ellAD>max(this.AdsorptionIsotherm_FT)) = 0;
    end

    function dP1D = interfacePotential(ell) 
        
        ellAD = min(this.AdsorptionIsotherm_FT)+ell;
        IP    = barychebevalMatrix(this.AdsorptionIsotherm_FT,ellAD);
        dP1D  = IP*IntPot;
        dP1D(ellAD<min(this.AdsorptionIsotherm_FT)) = 0;
     %   dP1D(ellAD>max(this.AdsorptionIsotherm_FT)) = 0;
    end

end