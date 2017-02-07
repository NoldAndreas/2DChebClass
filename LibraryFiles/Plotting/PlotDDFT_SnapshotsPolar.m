function PlotDDFT_SnapshotsPolar(input)

    %*****************************
    %Initialization of Data     
    v2struct(input);
    v2struct(data);
    
    [N1,N2,h1,h2,yrPlot,ytPlot,plotTimes] = LoadNumData(optsNum);    
        
    y1Plot     = Interp.pts1.*cos(Interp.pts2);
    y2Plot     = Interp.pts1.*sin(Interp.pts2);
    
    y1M    = reshape(y1Plot,Interp.Nplot2,Interp.Nplot1);
    y2M    = reshape(y2Plot,Interp.Nplot2,Interp.Nplot1);            
    
    xl = [(min(y1Plot)-0.5) (max(y1Plot)+0.5)];
    yl = [(min(y2Plot)-0.5) (max(y2Plot)+0.5)];

    
    ma3 = (~Ind.bound  & (Pts.y1_kv >=  min(yrPlot)) & (Pts.y1_kv <=  max(yrPlot)) & (Pts.y2_kv <=  max(ytPlot)) & (Pts.y2_kv >=  min(ytPlot)));
    maB = (Ind.bound   & (Pts.y1_kv >=  min(yrPlot)) & (Pts.y1_kv <=  max(yrPlot)) & (Pts.y2_kv <=  max(ytPlot)) & (Pts.y2_kv >=  min(ytPlot)));    
    
    rho_ic = rho_t(:,1);

    fl_norm = 0.1*max(max(max(abs(flux_t))))*[cos(max(ytPlot))  sin(max(ytPlot))]; 
    
    %Choose times for snapshots
    noPlots = 8;
    n       = length(plotTimes);
    d       = ceil(n/noPlots);
    ct      = 1:n;
    mark    = (mod(ct,d) == 1);
    
    cols = 4;
    rows = ceil(sum(mark)/cols);
    
    %**************************************
    %Initialization of figure and screen, and movie
    close all;
    figure        
    
    screen_size = get(0, 'ScreenSize');
    %f1 = figure(1);
    %set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    %set(f1, 'Position', [0 0 825 550]);%3000 2000] );
    
    set(gcf,'Color','white'); %Set background color    
    
    %*****************************    
    
    pl_j = 1;
    for j = 1:d:n
        
        x       = X_t(:,j);
        rho     = rho_t(:,j);
        t       = plotTimes(j);
        flux    = flux_t(:,j);
        
        
        subplot(rows,cols,pl_j);
               
        NormQuiverPlotPolar(Pts,flux,ma3,fl_norm,[0 0]); hold on;
        NormQuiverPlotPolar(Pts,flux,maB,fl_norm,[0 0],1,'k'); hold on;
        
       
        z       = real(Interp.InterPol*rho);                
    
        [C,h]   = contour(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1),[1 3 5 7 9],'LineWidth',1.0);    hold on;        
%        clabel(C,h,'FontSize',14);    

        
        pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1]);
        xlim(xl); ylim(yl);
        xlabel('y_1');        ylabel('y_2');
        title(['t = ',num2str(t)]);
        pl_j = pl_j +1;        
        
    end
    
        
    
end