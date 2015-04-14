function PlotDensitySlicesNormalInterface(this)

% Initialization         
    n      = 5;
    deltaZ = 5;
    [om,rho1D,params] = Compute1D(this,'LG');
    
    y2Max = this.optsNum.PlotAreaCart.y2Max;
    y1Min = this.optsNum.PlotAreaCart.y1Min;
    y1Max = this.optsNum.PlotAreaCart.y1Max;
        
    y1P = y1Min + (y1Max-y1Min)*(0.5:1:(n-0.5))/n;
    col = distinguishable_colors_NoRedBlueGreen();
    
    h_full  = this.hIII;
    dhIIIdx = this.y1_SpectralLine.Diff.Dy*h_full;
    %str = {'','','m','','c'};
        
    % Plotting
    figure('color','white','Position',[0 0 800 800]);        
    for i = 1:n
        y1              = y1P(i);
        IP              = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1P(i));
        h               = IP*h_full;        
        dhdx            = IP*dhIIIdx;
        [pts_y1{i},pts_y2{i}] = GetStartEndNormalPts(y1,h,dhdx,deltaZ);
        this.IDC.plotLine(pts_y1{i},0.5+pts_y2{i},this.GetRhoEq,struct('dist0',true,'plain',true,'CART',true,'color',col(i,:)));  hold on; %'plain',true
    end
    plot(this.IDC.Pts.y1*sin(this.IDC.alpha)+deltaZ-0.1,rho1D,'k--','linewidth',1.5);
    xlim([0 2]*deltaZ);
    
    
    figure('color','white','Position',[0 0 600 500]);    
    if(this.optsNum.PhysArea.alpha_deg ~= 90)
        PlotContourResults(this,{'hI','hII','hIII'}); hold on;        
    else
        PlotContourResults(this,{'hI','hII'}); hold on;        
    end
    for i = 1:n
        plot(pts_y1{i},pts_y2{i},':','color',col(i,:),'linewidth',1.5);
    end        
    
    function [pts_y1,pts_y2] = GetStartEndNormalPts(y1,h,dhdx,deltaZ)
        alpha  = atan(dhdx);
        
        d_lower = h/cos(alpha);
        d_int   = min(d_lower,deltaZ);
        
        pts_y1 = y1 + [-deltaZ,d_int]*sin(alpha);
        pts_y2 = h +  [deltaZ,-d_int]*cos(alpha);
    end
        
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