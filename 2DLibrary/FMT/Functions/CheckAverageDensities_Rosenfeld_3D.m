 function CheckAverageDensities_Rosenfeld_3D(this,IntMatrFex_2D,VisualOutput)           
 
       if(nargin == 2)
           VisualOutput = true;
       end
         
        IntMatrFex  = Get1DMatrices(IntMatrFex_2D,this);
        R = this.R(1);
        subPts.y2_kv              = (0.:0.01:3.5)';
        subPts.y1_kv              = inf*ones(size(subPts.y2_kv));           
        IP                        = this.AD.SubShapePts(subPts);
        Interp1D_AD.InterPol      = IP(:,this.AD.Pts.y1_kv == inf);
        Interp1D_AD.pts1          = subPts.y1_kv; 
        Interp1D_AD.pts2          = subPts.y2_kv;
        Interp1D_AD.ptsCart       = this.GetCartPts(Interp1D_AD.pts1,Interp1D_AD.pts2);

        subPtsAAD.y2_kv           = (this.y2Min:0.01:3.5)';
        subPtsAAD.y1_kv           = inf*ones(size(subPtsAAD.y2_kv));           
        IP                        = this.SubShapePts(subPtsAAD);
        Interp1D_AAD.InterPol     = IP(:,this.Pts.y1_kv == inf);
        Interp1D_AAD.pts1         = subPtsAAD.y1_kv; 
        Interp1D_AAD.pts2         = subPtsAAD.y2_kv;
        Interp1D_AAD.ptsCart      = this.GetCartPts(Interp1D_AAD.pts1,Interp1D_AAD.pts2);
        
        ADCartPts                 = this.AD.GetCartPts();
        CartPts                   = this.GetCartPts();
        
        markInf                   = (this.Pts.y1_kv== inf);
        markInfAD                 = (this.AD.Pts.y1_kv == inf);
 
        y0 = ones(this.N2,1);
        
        x = 0:0.01:5;
        
        if(VisualOutput)
            f1 = figure('name','Average densities');
            set(gcf,'Color','white');
            set(f1, 'Position', [0 0 500 1200]);
        end
        
        rows = 4;
        %*************************************************
        %1st Check: Average Densities with density = 1        
        if(VisualOutput)
            PlotRosenfeldFMT_AverageDensitiesInf(IntMatrFex(1),y0,rows,1); 
            subplot(rows,3,1); hold on; plot(x,check1n2(x),'m');
            subplot(rows,3,2); hold on; plot(x,check1n3(x),'m');
            subplot(rows,3,3); hold on; plot(x,check1n2z(x),'m');      
        end
        
        xADInt = ADCartPts.y2_kv(this.AD.Pts.y1_kv == inf);
        PrintErrorPos(IntMatrFex.AD.n2*y0-check1n2(xADInt),'FMT n2 for ones',xADInt);
        PrintErrorPos(IntMatrFex.AD.n3*y0-check1n3(xADInt),'FMT n3 for ones',xADInt);
        PrintErrorPos(IntMatrFex.AD.n2_v_2*y0-check1n2z(xADInt),'FMT n2_v_2 for ones',xADInt);
        %*************************************************
        
        %*************************************************
        %2nd Check: Average Densities with density = 2-erf(x)
        x = CartPts.y2_kv(markInf);    
        row = 2;
        if(VisualOutput)
            PlotRosenfeldFMT_AverageDensitiesInf(IntMatrFex(1),2-erf(x),rows,4);
            xC = 0:0.01:5;        
            subplot(rows,3,4); plot(xC,f2_h2(xC+R) - f2_h2(max(R,xC-R)),'m');%n2        
            subplot(rows,3,5); plot(xC,f3_h2(xC,xC+R) - f3_h2(xC,max(R,xC-R)),'m'); %n3       
            subplot(rows,3,6); plot(xC,f2z_h2(xC,xC+R) - f2z_h2(xC,max(R,xC-R)),'m');%n2z
        end
        
        PrintErrorPos(IntMatrFex.AD.n2*(2-erf(x))-f2_h2Diff(xADInt),'FMT n2 for 2-erf(x)',xADInt); %(f2_h2(xADInt+R) - f2_h2(max(R,xADInt-R)))
        PrintErrorPos(IntMatrFex.AD.n3*(2-erf(x))-f3_h2Diff(xADInt),'FMT n3 for 2-erf(x)',xADInt); %(f3_h2(xADInt,xADInt+R) - f3_h2(xADInt,max(R,xADInt-R)))
        PrintErrorPos(IntMatrFex.AD.n2_v_2*(2-erf(x))-f2z_h2Diff(xADInt),'FMT n2_v_2 for 2-erf(x)',xADInt); %(f2z_h2(xADInt,xADInt+R) - f2z_h2(xADInt,max(R,xADInt-R)))
        %*************************************************
        
        %*************************************************
        %3rd Check: Average Densities with density = 2-erf(x)
        z = ADCartPts.y2_kv(markInfAD);
        rhoAD = erf(z);
        row = 3;
        
        if(VisualOutput)
            PlotRosenfeldFMT_AADInf(IntMatrFex(1),rhoAD,rows,7);
            subplot(rows,3,6+1); plot(z,checkADn2(z),'m');            
            subplot(rows,3,6+2); plot(z,checkADn3(z),'m');
            subplot(rows,3,6+3); plot(z,checkADn2_v_2(z),'m');            
        end
        PrintErrorPos(IntMatrFex.AAD.n2*rhoAD-checkADn2(x),'FMT,AD n2 for ones',xADInt);
        PrintErrorPos(IntMatrFex.AAD.n3*rhoAD-checkADn3(x),'FMT,AD n3 for ones',xADInt);
        PrintErrorPos(IntMatrFex.AAD.n2_v_2*rhoAD-checkADn2_v_2(x),'FMT,AD n2_v_2 for ones',xADInt);
        
        %*************************************************
        
        
        %*********************************************
        z     = this.AD.Pts.y2_kv(this.AD.Pts.y1_kv == inf);    
        rhoAD = f2z_h2(z,z+R) - f2z_h2(z,max(R,z-R)); rhoAD(end) = 0;
        row   = 4;
        if(VisualOutput)
            PlotRosenfeldFMT_AADInf(IntMatrFex(1),rhoAD,rows,10);
        end
        %*********************************************
        
                
        %*******************************        
        zD = this.Pts.y2_kv(this.Pts.y1_kv == inf);    
        [FMT,phi] = Fex_FMTRosenfeld_3DFluid(ones(size(zD)),IntMatrFex,1,R);
        FMT_Bulk  = FexBulk_FMTRosenfeld_3DFluid(1,1);
        PrintErrorPos(FMT(end)-FMT_Bulk,'Bulk FMT Value');
%         if(abs(FMT(end)-FMT_Bulk) > 1e-3)            
%             cprintf('red','Error for computation of bulk FMT value > 10^(-3) \n');
%             figure('name','FMT');
%             do1Dplot_D(FMT);
%         else
%             cprintf('green','Error for computation of bulk FMT value < 10^(-3) \n');
%         end

%        cprintf('blue','Please check for errors in computation of average densities. \n');        
%        pause(2);       
%         close(f1);

    function y = f3_h2Diff(x)
        y          = (f3_h2(x,x+R) - f3_h2(x,max(R,x-R)));
        %y         = f3_h2(x,x+R) - f3_h2(max(R,x-R));
        y(x==Inf)  = 4/3*pi*R^3;
    end
    function y = f3_h2(x,xt)
        y = pi*((xt.^3.*erf(xt))/3 + xt.*(2*R^2 - 2*x.^2) + 2*x.*xt.^2 + (x.*erf(xt))/2 - (2*xt.^3)/3 + (xt.^2.*exp(-xt.^2))/(3*pi^(1/2)) - x.*xt.^2.*erf(xt) + (exp(-xt.^2).*(- R^2 + x.^2 + 1/3))/pi^(1/2) - xt.*erf(xt).*(R^2 - x.^2) - (x.*xt.*exp(-xt.^2))/pi^(1/2));
    end
 
    function y = f2_h2Diff(x)
        y          = f2_h2(x+R) - f2_h2(max(R,x-R));
        y(x==Inf)  = 2*pi*R;
    end
    function y = f2_h2(xt)
        y = pi*(2*xt - exp(-xt.^2)/pi^(1/2) - xt.*erf(xt));
    end

    function y = f2z_h2Diff(x)
        y          = f2z_h2(x,x+R) - f2z_h2(x,max(R,x-R));    
        y(x==Inf)  = 0;
    end
    function y = f2z_h2(x,xt)
        y = pi^(1/2)*xt.*exp(-xt.^2) - (pi*erf(xt))/2 - 2*pi^(1/2)*x.*exp(-xt.^2) - 2*pi*xt.^2 + pi*xt.^2.*erf(xt) + 4*pi*x.*xt - 2*pi*x.*xt.*erf(xt);
    end

    function y  = check1n2(x)
        y        = ones(size(x));
        y(x < 1) = x(x<2*R);
        y = y*2*pi*R;
    end
    function y  = check1n3(x)
        y        = 4*(R^3)*ones(size(x));
        y(x < 1) = (x(x<1)).^2.*(3*R-x(x<1));
        y        = y*pi/3;
    end
    function y  = check1n2z(x)
        y        = zeros(size(x));
        y(x < 1) = pi*(x(x<1)).*(x(x<1)-2*R);        
    end

     function y = checkADn2(z)
         y = -2*sqrt(pi)*R*(-z*sqrt(pi).*erf(z+R)-sqrt(pi)*erf(z+R)*R-exp(-(z+R).^2)-sqrt(pi)*z.*erf(-z+R)+sqrt(pi)*erf(-z+R)*R+exp(-(-z+R).^2));
         y(z==inf) = 4*pi*R^2;
     end
    function y = checkADn3(z) 
        y = (1/6)*sqrt(pi)*(-2*sqrt(pi)*erf(z+R).*z.^3.*exp(z.^2+2*R*z+R^2)+6*sqrt(pi)*erf(z+R).*z*R^2.*exp(z.^2+2*R*z+R^2)+4*sqrt(pi)*erf(z+R)*R^3.*exp(z.^2+2*R*z+R^2)-2*z.^2+2*R*z+4*R^2-2-3*sqrt(pi)*erf(z+R).*z.*exp(z.^2+2*R*z+R^2)-2*exp(z.^2+2*R*z+R^2).*erf(-z+R).*z.^3*sqrt(pi)+6*exp(z.^2+2*R*z+R^2).*erf(-z+R).*z*R^2*sqrt(pi)-4*exp(z.^2+2*R*z+R^2).*erf(-z+R)*R^3*sqrt(pi)+2*exp(4*R*z).*z.^2+2*exp(4*R*z)*R.*z-4*exp(4*R*z)*R^2+2*exp(4*R*z)-3*exp(z.^2+2*R*z+R^2).*z.*erf(-z+R)*sqrt(pi)).*exp(-z.^2-2*R*z-R^2);
        y(z>24) = 4/3*pi*R^3; %error < 10^-11
    end
    function y = checkADn2_v_2(z) 
        y = -(1/2)*sqrt(pi)*(-2*erf(z+R)*sqrt(pi).*exp(2*z.^2+2*R^2).*z.^2+2*erf(z+R)*sqrt(pi).*exp(2*z.^2+2*R^2)*R^2-2*exp((-z+R).^2).*z+2*exp((-z+R).^2)*R-erf(z+R)*sqrt(pi).*exp(2*z.^2+2*R^2)-2*erf(-z+R)*sqrt(pi).*exp(2*z.^2+2*R^2).*z.^2+2*erf(-z+R)*sqrt(pi).*exp(2*z.^2+2*R^2)*R^2+2*exp((z+R).^2).*z+2*exp((z+R).^2)*R-erf(-z+R)*sqrt(pi).*exp(2*z.^2+2*R^2)).*exp(-2*z.^2-2*R^2);
        y(z>15) = 0; %error: < 10^-13
    end

     function PlotRosenfeldFMT_AverageDensitiesInf(FMTMatrices,rho,rows,no)

        if(nargin == 2)
            f1 = figure('name','Average densities');
            set(gcf,'Color','white');
            set(f1, 'Position', [0 0 1300 400]);
            rows = 1;
            no   = 1;
        end
        
        nStruct = FMTMatrices.AD;            
        
        subplot(rows,3,no);            
        do1Dplot(nStruct.n2*rho); 
        hh = title('$n_2$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',20);
        pbaspect([1 1 1]);

        subplot(rows,3,no+1);
        do1Dplot(nStruct.n3*rho);  
        hh = title('$n_3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',20);
        pbaspect([1 1 1]);

        subplot(rows,3,no+2);
        do1Dplot(nStruct.n2_v_2*rho);  
        hh = title('${\bf n}_{2,y}$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',20);  
        pbaspect([1 1 1]);

    end

    function PlotRosenfeldFMT_AADInf(FMTMatrices,rho,rows,no)
        
        if(nargin == 2)
            f1 = figure('name','AAD');
            set(gcf,'Color','white');
            set(f1, 'Position', [0 0 1300 400]);
            row = 1;
        end
        
        nStruct = FMTMatrices.AAD;                
        
        subplot(rows,3,no); 
        do1Dplot_D(nStruct.n2*rho); 
        hh = title('$AAD_2$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',20);
        pbaspect([1 1 1]);
        
        subplot(rows,3,no+1); 
        do1Dplot_D(nStruct.n3*rho); title('n3');    
        hh = title('$AAD_3$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',20);
        pbaspect([1 1 1]);
        
        subplot(rows,3,no+2); 
        do1Dplot_D(nStruct.n2_v_2*rho); title('n2_2');    
        hh = title('$AAD,{\bf n}_{2,y}$'); set(hh,'Interpreter','Latex'); set(hh,'fontsize',20);
        pbaspect([1 1 1]);

    end

    function do1Dplot_D(val)
            plot(CartPts.y2_kv(markInf),val,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); hold on
            plot(Interp1D_AAD.ptsCart.y2_kv,Interp1D_AAD.InterPol*val,'linewidth',1.5);
            xlim([min(Interp1D_AAD.ptsCart.y2_kv) max(Interp1D_AAD.ptsCart.y2_kv)]);    
            h = xlabel('$y$'); set(h,'Interpreter','Latex'); set(h,'fontsize',25);        	
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);                        
    end
    function do1Dplot(val)
            plot(ADCartPts.y2_kv(this.AD.Pts.y1_kv == inf),val,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); hold on; %this.AD.Pts.y2_kv
            plot(Interp1D_AD.ptsCart.y2_kv,Interp1D_AD.InterPol*val,'linewidth',1.5);
            xlim([min(Interp1D_AD.ptsCart.y2_kv) max(Interp1D_AD.ptsCart.y2_kv)]);    
            h = xlabel('$y$');     set(h,'Interpreter','Latex'); set(h,'fontsize',25);
        	%h = ylabel('$\rho$');  set(h,'Interpreter','Latex'); set(h,'fontsize',25);                        
            %title(sel(i));    
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);                        
    end

end 