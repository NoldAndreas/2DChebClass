 function res = CheckAverageDensities_Rosenfeld_3D(this,IntMatrFex_2D,opts)           
  
       if((nargin > 2) && islogical(opts))
           VisualOutput = opts;
       else
           VisualOutput = true;
       end
       
       if(isa(this,'InfCapillary'))
           topB = true;                   
       else
           topB = false;
       end
       
       if(isa(this,'InfCapillary_FMT'))
           y2Min    = this.y2Min;
           checkTop = true;
       else
           y2Min    = 0.5;
           checkTop = false;
       end
         
        IntMatrFex  = Get1DMatrices(IntMatrFex_2D,this);
        R = this.R(1);
                
        if(topB)
            subPts.y2_kv          = (y2Min-R:0.01:this.y2Max+R)';
            subPtsAAD.y2_kv       = (y2Min:0.01:this.y2Max)';
        else
            subPts.y2_kv          = (0.:0.01:3.5)';
            subPtsAAD.y2_kv       = (this.y2Min:0.01:3.5)';        
        end
        subPts.y1_kv              = inf*ones(size(subPts.y2_kv));           
        IP                        = this.AD.SubShapePts(subPts);
        Interp1D_AD.InterPol      = IP(:,this.AD.Pts.y1_kv == inf);
        Interp1D_AD.pts1          = subPts.y1_kv; 
        Interp1D_AD.pts2          = subPts.y2_kv;
        Interp1D_AD.ptsCart       = this.GetCartPts(Interp1D_AD.pts1,Interp1D_AD.pts2);
        
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
        
        if(topB)
            x = (0:0.01:this.y2Max+R)';
        else            
            x = 0:0.01:5;
        end
        
        xADInt = ADCartPts.y2_kv(this.AD.Pts.y1_kv == inf);
        
        if(VisualOutput)
            f1 = figure('name','Average densities');
            set(gcf,'Color','white');
            set(f1, 'Position', [0 0 500 1200]);
        end
        
        rows = 4;
        %*************************************************
        %% 1st Check: Average Densities with density = 1        
        if(VisualOutput)
            PlotRosenfeldFMT_AverageDensitiesInf(IntMatrFex(1),y0,rows,1); 
            subplot(rows,3,1); hold on; plot(x,check1n2(x),'m');
            subplot(rows,3,2); hold on; plot(x,check1n3(x),'m');
            subplot(rows,3,3); hold on; plot(x,check1n2z(x),'m');      
        end
                
        res.error_n2_1   = PrintErrorPos(IntMatrFex.AD.n2*y0-check1n2(xADInt),'FMT n2 for ones',xADInt);
        res.error_n3_1   = PrintErrorPos(IntMatrFex.AD.n3*y0-check1n3(xADInt),'FMT n3 for ones',xADInt);
        res.error_n2v2_1 = PrintErrorPos(IntMatrFex.AD.n2_v_2*y0-check1n2z(xADInt),'FMT n2_v_2 for ones',xADInt);        
               
        %% 2nd Check: Average Densities with density = 2-erf(x)
        x       = CartPts.y2_kv(markInf);                
        
        if(checkTop)            
            xC_Cart = (this.y2Max-5):0.01:(this.y2Max+R); 
            x  = this.y2Max + R - x;    
            xC = this.y2Max + R - xC_Cart;                
        else
            xC_Cart = 0:0.01:5; 
            xC = xC_Cart;
        end
        row = 2;
        xADIntC = xADInt;
        if(topB && ~ checkTop)
            mark = (xADInt < 3);
        elseif(topB && checkTop)
            mark = (abs(this.y2Max + R - xADInt)  < 3);
            xADInt = this.y2Max + R - xADInt;
        else
            mark = true(size(xADInt));
        end
        
        rhoAD = 2-erf(x);
        
        if(VisualOutput)
            PlotRosenfeldFMT_AverageDensitiesInf(IntMatrFex,rhoAD,rows,4);                   
            subplot(rows,3,4); plot(xC_Cart,f2_h2(xC+R) - f2_h2(max(R,xC-R)),'m');%n2        
            subplot(rows,3,5); plot(xC_Cart,f3_h2(xC,xC+R) - f3_h2(xC,max(R,xC-R)),'m'); %n3       
            subplot(rows,3,6); plot(xC_Cart,f2z_h2(xC,xC+R) - f2z_h2(xC,max(R,xC-R)),'m');%n2z
        end
        
        err = IntMatrFex.AD.n2*(2-erf(x))-f2_h2Diff(xADInt);
        PrintErrorPos(err(mark),'FMT n2 for 2-erf(x)',xADIntC(mark)); 
        err = IntMatrFex.AD.n3*(2-erf(x))-f3_h2Diff(xADInt);
        PrintErrorPos(err(mark),'FMT n3 for 2-erf(x)',xADIntC(mark)); 
        err = IntMatrFex.AD.n2_v_2*(2-erf(x))-f2z_h2Diff(xADInt);
        PrintErrorPos(err(mark),'FMT n2_v_2 for 2-erf(x)',xADIntC(mark));
                                        
        %% 3rd Check: Average Densities with density = 2-erf(x)
        z_Cart    = ADCartPts.y2_kv(markInfAD);                                
        %xC_Cart = 0:0.01:5;
        if(checkTop)     
            z  = this.y2Max + R - z_Cart;
            %xC = this.y2Max + R - xC_Cart;                
        else
            z = z_Cart;
        end
        
        rhoAD = erf(z);           
        row   = 3;
        
        if(VisualOutput)
            PlotRosenfeldFMT_AADInf(IntMatrFex(1),rhoAD,rows,7);
            subplot(rows,3,6+1); plot(z_Cart,checkADn2(z),'m');            
            subplot(rows,3,6+2); plot(z_Cart,checkADn3(z),'m');
            subplot(rows,3,6+3); plot(z_Cart,checkADn2_v_2(z),'m');            
        end
        res.error_n2AD_erf = PrintErrorPos(IntMatrFex.AAD.n2*rhoAD-checkADn2(x),'FMT,AD n2 for erf',x);
        res.error_n3AD_erf = PrintErrorPos(IntMatrFex.AAD.n3*rhoAD-checkADn3(x),'FMT,AD n3 for erf',x);
        res.error_n2v2AD_erf = PrintErrorPos(IntMatrFex.AAD.n2_v_2*rhoAD-checkADn2_v_2(x),'FMT,AD n2_v_2 for erf',x);
        
        %*************************************************        
        %% 4th Check: Average Densities with density = 1 for y > 1.
                
       % z     = this.AD.Pts.y2_kv(this.AD.Pts.y1_kv == inf);    
        z     = ADCartPts.y2_kv(this.AD.Pts.y1_kv == inf);    
        if(checkTop)     
            z = this.y2Max + R - z;
        end
        rhoAD = ones(size(z));
        rhoAD(this.AD.mark_id_2(:,1)) = 0;%f2z_h2(z,z+R) - f2z_h2(z,max(R,z-R)); rhoAD(end) = 0;        
        if(checkTop)     
            rhoAD(this.AD.mark_id_2(:,3)) = 0;%f2z_h2(z,z+R) - f2z_h2(z,max(R,z-R)); rhoAD(end) = 0;
        end
        row   = 4;
        %xADInt =  this.Pts.y2_kv(this.Pts.y1_kv == inf) - 0.5;
        xADInt =  CartPts.y2_kv(this.Pts.y1_kv == inf);
        if(VisualOutput)
            PlotRosenfeldFMT_AADInf(IntMatrFex(1),rhoAD,rows,10);
            subplot(rows,3,9+1); plot(xADInt,check1n2_AD(xADInt),'m','linewidth',1.5);
            subplot(rows,3,9+2); plot(xADInt,check1n3_AD(xADInt),'m','linewidth',1.5);
            subplot(rows,3,9+3); plot(xADInt,check1n2z_AD(xADInt),'m','linewidth',1.5);            
        end        
        res.error_n2AD_ones   = PrintErrorPos(IntMatrFex.AAD.n2*rhoAD-check1n2_AD(xADInt),'FMT, AD n2 for ones',xADInt);
        res.error_n3AD_ones   = PrintErrorPos(IntMatrFex.AAD.n3*rhoAD-check1n3_AD(xADInt),'FMT, AD n3 for ones',xADInt);
        res.error_n2v2AD_ones = PrintErrorPos(IntMatrFex.AAD.n2_v_2*rhoAD-check1n2z_AD(xADInt),'FMT, AD n2_v_2 for ones',xADInt);        
        
        %*********************************************
        
                
        %*******************************       
        if(~topB)
            zD = this.Pts.y2_kv(this.Pts.y1_kv == inf);    
            [FMT,phi] = Fex_FMTRosenfeld_3DFluid(ones(size(zD)),IntMatrFex,1,R);
            FMT_Bulk  = FexBulk_FMTRosenfeld_3DFluid(1,1);
            PrintErrorPos(FMT(end)-FMT_Bulk,'Bulk FMT Value');
        end


    function y = f3_h2Diff(x)
        y          = (f3_h2(x,x+R) - f3_h2(x,max(R,x-R)));
        y(x > 20)  = 4/3*pi*R^3;
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
        y(x>20)  = 0;                
    end
    function y = f2z_h2(x,xt)
        y = pi^(1/2)*xt.*exp(-xt.^2) - (pi*erf(xt))/2 - 2*pi^(1/2)*x.*exp(-xt.^2) - 2*pi*xt.^2 + pi*xt.^2.*erf(xt) + 4*pi*x.*xt - 2*pi*x.*xt.*erf(xt);
        if(topB && checkTop)
            y = -y;
        end
    end

    function y  = check1n2_AD(x)
        y        = ones(size(x));
        
        markB    = (x < (y2Min + 2*R));
        y(markB) =  x(markB) - (y2Min);
        
        if(topB)
            markT    = ((this.y2Max -x) < 1);
            y(markT) = this.y2Max - x(markT);
        end
        y = y*2*pi*R;
    end
    function y  = check1n2(x)
        y        = ones(size(x));
        
        markB    = (x < (y2Min+R));
        y(markB) =  x(markB) - (y2Min - R);                
        
        if(topB)
            markT    = ((this.y2Max+R -x) < 1);
            y(markT) = this.y2Max + R - x(markT);
        end
        y = y*2*pi*R;
    end

    function y  = check1n3_AD(x)
        y        = 4*(R^3)*ones(size(x));
               
        markB    = (x < (y2Min+2*R));
        xR       = x(markB) - (y2Min);
        y(markB) = xR.^2.*(3*R-xR);            
        %y(x < 1) = (x(x<1)).^2.*(3*R-x(x<1));
        
        if(topB)
            markT    = ((this.y2Max -x) < 1);                        
            xR       = this.y2Max - x(markT);            
            y(markT) = xR.^2.*(3*R-xR);                        
        end
        
        y        = y*pi/3;
    end
    function y  = check1n3(x)
        y        = 4*(R^3)*ones(size(x));
               
        markB    = (x < (y2Min+R));
        xR       = x(markB) - (y2Min - R);
        y(markB) = xR.^2.*(3*R-xR);            
        %y(x < 1) = (x(x<1)).^2.*(3*R-x(x<1));
        
        if(topB)
            markT    = ((this.y2Max+R -x) < 1);
            xR       = this.y2Max + R - x(markT);            
            y(markT) = xR.^2.*(3*R-xR);                        
        end
        
        y        = y*pi/3;
    end

    function y  = check1n2z_AD(x)
        y        = zeros(size(x));
        
        markB    = (x < (y2Min+2*R));
        xR       = x(markB) - (y2Min);
                
        y(markB) = pi*xR.*(xR-2*R);                    
        %y(x < 1) = pi*(x(x<1)).*(x(x<1)-2*R);        
        
        if(topB)
            markT    = ((this.y2Max -x) < 1);
            xR       = this.y2Max - x(markT);           
            y(markT) = -pi*xR.*(xR-2*R);                    
        end
        
    end
    function y  = check1n2z(x)
        y        = zeros(size(x));        
        markB    = (x < (y2Min+R));
        xR       = x(markB) - (y2Min - R);
        y(markB) = pi*xR.*(xR-2*R);                    
        
        if(topB)
            markT    = ((this.y2Max+R -x) < 1);
            xR       = this.y2Max + R - x(markT);           
            y(markT) = -pi*xR.*(xR-2*R);                    
        end
        
    end

    function y = checkADn2(z)        
         %y = -2*sqrt(pi)*R*(-z*sqrt(pi).*erf(z+R)-sqrt(pi)*erf(z+R)*R-exp(-(z+R).^2)-sqrt(pi)*z.*erf(-z+R)+sqrt(pi)*erf(-z+R)*R+exp(-(-z+R).^2));                  
         y = -2*pi*R*(-z.*erf(z+R)-erf(z+R)*R-exp(-(z+R).^2)/sqrt(pi)-(z-R).*erf(-z+R)+exp(-(-z+R).^2)/sqrt(pi));                           
         y(z==inf) = 4*pi*R^2;                           
     end
    function y = checkADn3(z) 
        y = (1/6)*sqrt(pi)*(-2*sqrt(pi)*erf(z+R).*z.^3.*exp(z.^2+2*R*z+R^2)+6*sqrt(pi)*erf(z+R).*z*R^2.*exp(z.^2+2*R*z+R^2)+4*sqrt(pi)*erf(z+R)*R^3.*exp(z.^2+2*R*z+R^2)-2*z.^2+2*R*z+4*R^2-2-3*sqrt(pi)*erf(z+R).*z.*exp(z.^2+2*R*z+R^2)-2*exp(z.^2+2*R*z+R^2).*erf(-z+R).*z.^3*sqrt(pi)+6*exp(z.^2+2*R*z+R^2).*erf(-z+R).*z*R^2*sqrt(pi)-4*exp(z.^2+2*R*z+R^2).*erf(-z+R)*R^3*sqrt(pi)+2*exp(4*R*z).*z.^2+2*exp(4*R*z)*R.*z-4*exp(4*R*z)*R^2+2*exp(4*R*z)-3*exp(z.^2+2*R*z+R^2).*z.*erf(-z+R)*sqrt(pi)).*exp(-z.^2-2*R*z-R^2);
        y(z>24) = 4/3*pi*R^3; %error < 10^-11
    end
    function y = checkADn2_v_2(z) 
        y = -(1/2)*sqrt(pi)*(-2*erf(z+R)*sqrt(pi).*exp(2*z.^2+2*R^2).*z.^2+2*erf(z+R)*sqrt(pi).*exp(2*z.^2+2*R^2)*R^2-2*exp((-z+R).^2).*z+2*exp((-z+R).^2)*R-erf(z+R)*sqrt(pi).*exp(2*z.^2+2*R^2)-2*erf(-z+R)*sqrt(pi).*exp(2*z.^2+2*R^2).*z.^2+2*erf(-z+R)*sqrt(pi).*exp(2*z.^2+2*R^2)*R^2+2*exp((z+R).^2).*z+2*exp((z+R).^2)*R-erf(-z+R)*sqrt(pi).*exp(2*z.^2+2*R^2)).*exp(-2*z.^2-2*R^2);
        y(z>15) = 0; %error: < 10^-13
        
        if(topB && checkTop)
            y = -y;
        end
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
            set(gca,'fontsize',20);                        
            set(gca,'linewidth',1.5);                        
    end

end 