function PlotDensitySlicesNormalInterface(this,y1P)

% Initialization             
    deltaZ = 5;
    [~,rho1D_LV] = Compute1D(this,'LG');
    if(this.optsNum.PhysArea.alpha_deg > 90)
       rho1D_LV = flipud(rho1D_LV);
    end
    
    
    [~,rho1D_WL] = Compute1D(this,'WL');
    rho1D_WL_full = kron(ones(this.IDC.Pts.N1,1),rho1D_WL);
    
    
    %y2Max = this.optsNum.PlotAreaCart.y2Max;
    if(nargin < 2)
        n      = 5;
        y1Min = this.optsNum.PlotAreaCart.y1Min;
        y1Max = this.optsNum.PlotAreaCart.y1Max;
        y1P = y1Min + (y1Max-y1Min)*(0.5:1:(n-0.5))/n;    
    end
    
    rhoLiq_sat = this.optsPhys.rhoLiq_sat;
        
    
    col = distinguishable_colors_NoRedBlueGreen();
    
    h_full  = this.hIII;
    dhIIIdx = this.y1_SpectralLine.Diff.Dy*h_full;
    %str = {'','','m','','c'};
                
    L2_AdIso = 4;%input('Parameter L_2 used for computation of adsorption isotherm: ');
    for i = 1:length(y1P)
        
        if( abs(90 - this.alpha_YCA*180/pi) > 5)
            y1              = y1P(i);
            IP              = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1P(i));
            h               = IP*h_full;        
            dhdx            = IP*dhIIIdx;
            [pts_y1{i},pts_y2{i}] = GetStartEndNormalPts(y1,h,dhdx,deltaZ);
        else
            pts_y1{i} = deltaZ*[-1,1];
            pts_y2{i} = [1,1]*this.optsNum.PlotAreaCart.y2Max*(0.5+(i-1))/n;
        end
        [f_p{i},pts{i}]       = this.IDC.plotLine(pts_y1{i},0.5+pts_y2{i},this.GetRhoEq,struct('dist0',true,'plain',true,'CART',true,'color',col(i,:))); %'plain',true        
        
        %rho_Ana = rholv(z)*rho_wl(y_2)/rho_liq
        IP_WL      = this.IDC.SubShapePtsCart(pts{i});                
        rho1D_wl_i = IP_WL*rho1D_WL_full;
                
        IP_LV      = barychebevalMatrix(this.IDC.Pts.x1,this.IDC.CompSpace1((pts{i}.z-deltaZ)/sin(this.IDC.alpha)));
        rho1D_lv_i = IP_LV*rho1D_LV;
        rho_Ana{i} = (rho1D_lv_i.*rho1D_wl_i)/rhoLiq_sat;
        %plot(pts{i}.z-deltaZ,rho_Ana{i},'--','color',col(i,:),'linewidth',1.5); hold on;       
        
        %rho_Ana2 = rho_dp(y_2)
        if(~isempty(this.AdsorptionIsotherm) && ( abs(90 - this.alpha_YCA*180/pi) > 5))            
            th       = atan(dhdx);
            ell      = this.y1_SpectralLine.InterpolationMatrix_Pointwise(y1P(i))*this.hIII/cos(th);        
            [rho,mu] = GetPointAdsorptionIsotherm(this,ell);        
            y2_AdIso = this.AdsorptionIsotherm.Pts.y2;
            x2_AdIso = this.AdsorptionIsotherm.Pts.x2;

            x2       = InvQuotientMap((pts{i}.y2_kv-0.5)/cos(th)+0.5,L2_AdIso,0.5,inf);
            IP_AdIso = barychebevalMatrix(x2_AdIso,x2);
            rho_Ana2{i} = IP_AdIso*rho';
        end
        %plot(y2_AdIso,rho); hold on;
        %plot(this.IDC.Pts.y1*sin(this.IDC.alpha)+ell+0.5,flipud(rho1D_LV),'m')
        %plot(this.AdsorptionIsotherm.Pts.y2-0.5,rho,'--','color',col(i,:),'linewidth',1.5); hold on;%%%%
    end
    
    close all
    
    figure('color','white','Position',[0 0 600 600]);
    
    for i = 1:length(y1P)
        plot(pts{i}.z-deltaZ,f_p{i},'color',col(i,:)); hold on;
        if(~isempty(this.AdsorptionIsotherm))
            %plot(pts{i}.z-deltaZ,rho_Ana1{i},'--','color',col(i,:),'linewidth',1.5); hold on;       
            %plot(pts{i}.z-deltaZ,rho_Ana2{i},'--','color',col(i,:),'linewidth',1.5); hold on;       
        end
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
          
    SaveCurrentFigure(this,'DensityNormalInterface');          
        
    figure('color','white','Position',[0 0 600 600]);
    if(this.optsNum.PhysArea.alpha_deg ~= 90)
        PlotContourResults(this,{'hI','hII','hIII'}); hold on;        
    else
        PlotContourResults(this,{'hI','hII'}); hold on;        
    end
    for i = 1:length(y1P)
        plot(pts_y1{i},pts_y2{i},':','color',col(i,:),'linewidth',1.3);
    end        
    
    SaveCurrentFigure(this,'DensitySlicesNormal_Contour');
    
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