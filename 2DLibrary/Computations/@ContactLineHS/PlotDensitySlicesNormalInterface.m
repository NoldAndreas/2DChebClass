function PlotDensitySlicesNormalInterface(this)

% Initialization         
    n      = 5;
    deltaZ = 5;
    [~,rho1D_LV] = Compute1D(this,'LG');
    if(this.optsNum.PhysArea.alpha_deg > 90)
       rho1D_LV = flipud(rho1D_LV);
    end
    
    
    [~,rho1D_WL] = Compute1D(this,'WL');
    rho1D_WL_full = kron(ones(this.IDC.Pts.N1,1),rho1D_WL);
    
    
    y2Max = this.optsNum.PlotAreaCart.y2Max;
    y1Min = this.optsNum.PlotAreaCart.y1Min;
    y1Max = this.optsNum.PlotAreaCart.y1Max;
    
    rhoLiq_sat = this.optsPhys.rhoLiq_sat;
        
    y1P = y1Min + (y1Max-y1Min)*(0.5:1:(n-0.5))/n;
    col = distinguishable_colors_NoRedBlueGreen();
    
    h_full  = this.hIII;
    dhIIIdx = this.y1_SpectralLine.Diff.Dy*h_full;
    %str = {'','','m','','c'};
                
    for i = 1:n
        y1              = y1P(i);
        IP              = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1P(i));
        h               = IP*h_full;        
        dhdx            = IP*dhIIIdx;
        [pts_y1{i},pts_y2{i}] = GetStartEndNormalPts(y1,h,dhdx,deltaZ);
        [f_p{i},pts{i}] = this.IDC.plotLine(pts_y1{i},0.5+pts_y2{i},this.GetRhoEq,struct('dist0',true,'plain',true,'CART',true,'color',col(i,:))); %'plain',true        
        
        IP_WL      = this.IDC.SubShapePtsCart(pts{i});                
        rho1D_wl_i = IP_WL*rho1D_WL_full;
                
        IP_LV      = barychebevalMatrix(this.IDC.Pts.x1,this.IDC.CompSpace1((pts{i}.z-deltaZ)/sin(this.IDC.alpha)));
        rho1D_lv_i = IP_LV*rho1D_LV;
        rho_Ana{i} = (rho1D_lv_i.*rho1D_wl_i)/rhoLiq_sat;
        plot(pts{i}.z-deltaZ,rho_Ana{i},'--','color',col(i,:),'linewidth',1.5); hold on;       
    end
    
    close all
    
    figure('color','white','Position',[0 0 1700 600]);
    subplot(1,2,1);
    for i = 1:n
        plot(pts{i}.z-deltaZ,f_p{i},'color',col(i,:),'linewidth',1.5); hold on;
        plot(pts{i}.z-deltaZ,rho_Ana{i},'--','color',col(i,:),'linewidth',1.5); hold on;       
    end    
    dz = (-deltaZ:0.05:deltaZ)';
    IP_LV      = barychebevalMatrix(this.IDC.Pts.x1,this.IDC.CompSpace1(dz/sin(this.IDC.alpha)));
	rho1D_lv_i = IP_LV*rho1D_LV;
	plot(dz,rho1D_lv_i,'k--','linewidth',1.5); hold on;        
    %plot(this.IDC.Pts.y1*sin(this.IDC.alpha),rho1D_LV,'k--','linewidth',1.5);
    
    
    
    
    xlim([-1 1]*deltaZ);
    ylim([0 1.4]);
    set(gca,'fontsize',20);
    set(gca,'linewidth',1.5);
    xlabel('$z$','Interpreter','Latex','fontsize',25);
    ylabel('$n \sigma^3$','Interpreter','Latex','fontsize',25);
    pbaspect([(this.optsNum.PlotAreaCart.y1Max-this.optsNum.PlotAreaCart.y1Min) ....
              (this.optsNum.PlotAreaCart.y2Max-this.optsNum.PlotAreaCart.y2Min) 1]);
        
    subplot(1,2,2);
    if(this.optsNum.PhysArea.alpha_deg ~= 90)
        PlotContourResults(this,{'hI','hII','hIII'}); hold on;        
    else
        PlotContourResults(this,{'hI','hII'}); hold on;        
    end
    for i = 1:n
        plot(pts_y1{i},pts_y2{i},':','color',col(i,:),'linewidth',1.5);
    end        
    
    SaveCurrentFigure(this,'DensityNormalInterface');
    
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
% 
%    function [DeltaY1_II,DeltaY1_III] = ComputeDeltaFit(rhoLV)            
%         
%         dP1D        = GetDisjoiningPressure_I_ell(this,this.hI);
%         [min_I,i_I] = min(dP1D);
%         [min_III,i_III] = min(GetDisjoiningPressure_III(this));
%         DeltaY1_III = this.y1_SpectralLine.Pts.y(i_III) - this.y1_I(i_I);
% 
%         hS           = (this.optsPhys.rhoGas_sat + this.optsPhys.rhoLiq_sat)/2;        
%         fsolveOpts   = optimset('Display','off');            
%         
%      
% 
%         f            = rhoLV;
%         [~,j]        = max(this.hII);
%         [DeltaY1_II,~,exitflag] = fsolve(@fX,fsolveOpts);            
%         if(exitflag < 1)
%             cprintf('*r','ComputeDeltaFit: Fitting hII vs hI: no solution found');
%         end
% 
%         function z = fX(y1)                
%             IP_h = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1);
%             z    = IP_h*f-hS;
%         end    
% 	end

end