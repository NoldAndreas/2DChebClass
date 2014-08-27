function PlotDynamicContactAngles

    plotLegend = {};
    noPlot = 1;

    figure('color','white','Position',[0 0 800 800]);
        
    LoadDataFluidData('Stroem/Stroem_ParaffinOil_untreatedPS','s','r');   hold on;  
    %LoadDataFluidData('Stroem/Stroem_ParaffinOil_PTFE','d','r'); %(no zero static contact angle)
    LoadDataFluidData('Stroem/Stroem_SiliconeOil_I_untreatedPS','d','r');    
    LoadDataFluidData('Stroem/Stroem_SiliconeOil_II_untreatedPS','o','r');    
    LoadDataFluidData('Stroem/Stroem_SiliconeOil_II_oxidizedPS','v','r');    
    
 
    LoadHoffmannData('BlackCircles','s','k');
    LoadHoffmannData('Triangles','d','k');
    LoadHoffmannData('Crosses','+','k');
    LoadHoffmannData('Hexagons','o','w');
    LoadHoffmannData('Squares','s','w');            
    
    LoadHoffmannData('Fit','-',[]);    
       
    h_legend = legend(plotLegend,'Location','Northwest');          
    set(h_legend,'FontSize',8);
    
    thetaM = (0:0.01:pi)';
    G_thM  = GHR(thetaM);
    
    plot(G_thM/log(10^4),thetaM*180/pi,'k--','linewidth',2);
    
    xlim([1e-5 100]);
    xlabel('$Ca$','Interpreter','Latex','fontsize',20);
    ylabel('$\theta_m [^\circ]$','Interpreter','Latex','fontsize',20);
    
    set(gca,'fontsize',20);
    
    print2eps('Hoffmann_Data',gcf);  
	saveas(gcf,'Hoffmann_Data.fig');    
        
   
    function data =  LoadHoffmannData(name,symbol,color)
        fid = fopen(['D://Data/ExperimentalContactAngle/Hoffmann/Hoffmann_',name,'.txt']);
        y = textscan(fid,'%[^\n]',1,'headerlines',4); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f'); %[T, rhoG, rhoL]
        fclose(fid); 
        
        data.legend = char(y{1});
        data.Feq = x{1}(1);
        data.Ca  = x{1}(2:end) - data.Feq;
        data.theta = x{2}(2:end);
        
        if(isempty(color))
            semilogx(data.Ca,data.theta,symbol,'linewidth',2); hold on;
        else
            semilogx(data.Ca,data.theta,symbol,'MarkerFaceColor',color,'MarkerSize',10); hold on;
        end
        
        plotLegend{noPlot} = data.legend;
        noPlot             = noPlot + 1;
    end

    function data =  LoadDataFluidData(name,symbol,color)
        fid = fopen(['D://Data/ExperimentalContactAngle/',name,'.txt']);
        y = textscan(fid,'%[^\n]',1,'headerlines',4); %[T, rhoG, rhoL]        
        x = textscan(fid,'%f %f'); %[T, rhoG, rhoL]
        fclose(fid); 
        
        data.legend = char(y{1});
        data.gamma = x{1}(1);
        data.eta = x{1}(2);
        data.rho = x{1}(3);
        
        data.theta = x{1}(4:end);
        data.Ca    = x{2}(4:end)*data.eta/data.gamma*10^-3;
        
        if(isempty(color))
            semilogx(data.Ca,data.theta,symbol,'linewidth',2); hold on;
        else
            semilogx(data.Ca,data.theta,symbol,'MarkerFaceColor',color,'MarkerSize',10); hold on;
        end
        
        plotLegend{noPlot} = data.legend;
        noPlot             = noPlot + 1;
    end

    function z = GHR(t)
        z = 1i*pi^2/24 - t/2.*log(1+exp(1i*t)) ...
            + 1i/2*(  dilog(1+exp(1i*t)) + ...
                      dilog(exp(1i*t)) ) - sin(t)/2;
        disp(['Max imaginary part: ',num2str(max(imag(z)))]);
        z = real(z);
    end
end