function ConvergenceStudy(data)

    vext     = @Vext10; %Vext7
    simul    = @SimulationPolarInfinity;%Disk; %SimulationBox
    L1       = 2;

    n1 = 5:2:50;
    n2 = 4:2:20;
    
    if(nargin == 0)    
        for i=1:length(n1)
            for j=1:length(n2)                
                data(i,j) = simul(n1(i),n2(j),L1,1,vext);
            end 
        end
        SaveToFile('ConvergenceStudy',data,getResultsPath());
    end       
    
    n1D        = reshape([data.N1],length(n1),length(n2));
    n2D        = reshape([data.N2],length(n1),length(n2));
    dy1D       = reshape([data.dy1],length(n1),length(n2));
    ddy1D      = reshape([data.ddy1],length(n1),length(n2));
    dy2D       = reshape([data.dy2],length(n1),length(n2));
    ddy2D      = reshape([data.ddy2],length(n1),length(n2));    
    InterPolD  = reshape([data.InterPol],length(n1),length(n2));
    Int        = reshape([data.Int],length(n1),length(n2));        
%    convD       = reshape([data.Conv],length(n2),length(n1));    
    
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
    contourf(n1D,n2D,log(Int)/log(10)); xlabel('N_1'); ylabel('N_2'); 
    title('log_{10}(Error of Integration)');        
    colorbar;
    
    subplot(2,2,2); 
    
    contourf(n1D,n2D,log(InterPolD)/log(10)); xlabel('N_1'); ylabel('N_2'); 
    title('log_{10}(Error of Interpolation)'); 
    colorbar;
    
    subplot(2,2,3);
    semilogy(n1D(:),dy1D(:),'ok'); hold on;
    semilogy(n1D(:),ddy1D(:),'om'); 
    xlabel('N_1');
    legend('$\frac{df}{dy_1}$','$\frac{d^2f}{dy_1^2}$');
    h = legend;    set(h, 'interpreter', 'latex');
    LEG = findobj(h,'type','text');    set(LEG,'FontSize',15);
    %pbaspect([1 1 1]);
    %title('Numerical Errors');
    
    subplot(2,2,4);
    semilogy(n2D(:),dy2D(:),'og'); hold on;
    semilogy(n2D(:),ddy2D(:),'ob'); 
     xlabel('N_2');
    legend('$\frac{df}{dy_2}$','$\frac{d^2f}{dy_2^2}$')
    h = legend;    set(h, 'interpreter', 'latex');
    LEG = findobj(h,'type','text');    set(LEG,'FontSize',15);
    %pbaspect([1 1 1]);
    %title('Numerical Errors');       
   
    %title('$\frac{A-A(-1)}{Y}$','interpreter','latex')
    
    
end