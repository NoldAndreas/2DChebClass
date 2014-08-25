function PlotDynamicContactAngles

    figure('color','white','Position',[0 0 800 800]);
 
    LoadHoffmannData('BlackCircles','s','k');
    LoadHoffmannData('Triangles','d','k');
    LoadHoffmannData('Crosses','+','k');
    LoadHoffmannData('Hexagons','o','w');
    LoadHoffmannData('Squares','s','w');    
    
    LoadHoffmannData('Fit','-',[]);    
        
    thetaM = (0:0.01:pi)';
    G_thM  = GHR(thetaM);
    
    plot(G_thM/log(10^4),thetaM*180/pi,'k--','linewidth',2);
    
    xlim([1e-5 100]);
    xlabel('$Ca + F(\theta_{eq})$','Interpreter','Latex','fontsize',20);
    ylabel('$\theta_m [deg]$','Interpreter','Latex','fontsize',20);
    
    set(gca,'fontsize',20);
    
    print2eps('Hoffmann_Data',gcf);  
	saveas(gcf,'Hoffmann_Data.fig');    
        
   
    function x =  LoadHoffmannData(name,symbol,color)
        fid = fopen(['D://Data/ExperimentalContactAngle/Hoffmann_',name,'.txt']);
        x = textscan(fid,'%f %f','headerlines',4); %[T, rhoG, rhoL]
        fclose(fid); 
        
        if(isempty(color))
            semilogx(x{1},x{2},symbol,'linewidth',2); hold on;
        else
            semilogx(x{1},x{2},symbol,'MarkerFaceColor',color,'MarkerSize',10); hold on;
        end
    end

    function z = GHR(t)
        z = 1i*pi^2/24 - t/2.*log(1+exp(1i*t)) ...
            + 1i/2*(  dilog(1+exp(1i*t)) + ...
                      dilog(exp(1i*t)) ) - sin(t)/2;
        disp(['Max imaginary part: ',num2str(max(imag(z)))]);
        z = real(z);
    end
end