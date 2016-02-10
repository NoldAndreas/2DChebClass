function PlotInterfaceFittingQuality(this,nts)

    if(nargin < 2)
        y2Lim = [5 7.5];
    end
                
    y2_Interface= (this.optsNum.PlotAreaCart.y2Min:0.1:this.optsNum.PlotAreaCart.y2Max)';
    
    c       = 1/tan(this.optsNum.PhysArea.alpha_deg*pi/180);
    shapeSL = struct('yMin',this.optsNum.PlotAreaCart.y1Min,...
                         'yMax',this.optsNum.PlotAreaCart.y1Max,...
                         'N',40);                       
	this.y1_SpectralLine = SpectralLine(shapeSL);
    this.y1_SpectralLine.ComputeAll();
    y1 = this.y1_SpectralLine.Pts.y;
                  
    figure('color','white');
    if(nts(1) ~= 1)        
        nts = [1,nts];
    end
    
    for nt = nts            
        ExtrapolateInterface(this,this.dynamicsResult.rho_t(:,1,nt),...
                                   y2_Interface,...
                                   this.dynamicsResult.ak(:,nt),...
                                   [],{'analyse'},...
                                   this.dynamicsResult.t(nt));                                                                   
        hIII = Compute_hIII(this,this.dynamicsResult.rho_t(:,1,nt));
        
        mark = (hIII < 15);
        hIII = hIII(mark);
        y1P  = y1(mark);
        
        if(nt == 1);
            [~,calInd] = min(abs(hIII-10));            
            y10 = y1P(calInd) - c*hIII(calInd);            
        end        
        plot(hIII,y1P-c*hIII - y10,'b'); hold on;
        
        ax = gca;                
        set(ax,'YScale','log');
    end       
    pbaspect([1 1 1]);    
    SaveCurrentFigure(this,'InterfaceFitting');

end