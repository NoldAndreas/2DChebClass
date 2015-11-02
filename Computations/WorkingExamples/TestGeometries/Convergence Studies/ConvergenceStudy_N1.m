function ConvergenceStudy_N1(data)

    vext     = @Vext8; 
    simul    = @SimulationInfiniteCapillary; %@SimulationBox;%@SimulationPolarInfinity;
    L1       = 1;

    n1 = 10:4:40;
    n2 = 25;
    
    if(nargin == 0)    
        for i=1:length(n1)
            data(i) = simul(n1(i),n2,L1,1,vext);
        end
        SaveToFile('ConvergenceStudy_N1',data,getResultsPath());
    end       
    
    %************************************************************
    %*********************Plotting*******************************    
    %************************************************************
    
    close all;
    figure
    %screen_size = get(0, 'ScreenSize');
    f1 = figure(1);
    set(f1, 'Position', [0 0 1000 750 ] );    
    
    set(gcf,'Color','white'); %Set background color    
            
    subplot(2,2,1);     
    plot([data.N1],log([data.Int])/log(10)); xlabel('N_1'); ylabel('N_2'); 
    title('log_{10}(Error of Integration)');        
    
    subplot(2,2,2);     
    plot([data.N1],log([data.InterPol])/log(10)); xlabel('N_1'); ylabel('N_2'); 
    title('log_{10}(Error of Interpolation)'); 
    
    subplot(2,2,3);
    semilogy([data.N1],[data.dy1],'ok'); hold on;
    semilogy([data.N1],[data.ddy1],'om'); 
    xlabel('N_1');
    legend('$\frac{df}{dy_1}$','$\frac{d^2f}{dy_1^2}$');
    h = legend;    set(h, 'interpreter', 'latex');
    LEG = findobj(h,'type','text');    set(LEG,'FontSize',15);
    %pbaspect([1 1 1]);
    %title('Numerical Errors');
    
    
end