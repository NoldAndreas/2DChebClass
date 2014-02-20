function PlotDDFT_Snapshots_2Species(input)

    %*****************************
    %Initialization of Data     
    v2struct(input);
    v2struct(data);
    
    [N1,N2,h1,h2,y1Plot,y2Plot,plotTimes] = LoadNumData(optsNum);    
    
    xl = [(min(y1Plot)-0.5) (max(y1Plot)+0.5)];
    yl = [(min(y2Plot)-0.5) (max(y2Plot)+0.5)];

    
    ma3 = (~Ind.bound  & (Pts.y1_kv >=  min(y1Plot)) & (Pts.y1_kv <=  max(y1Plot)) & (Pts.y2_kv <=  max(y2Plot)) & (Pts.y2_kv >=  min(y2Plot)));
    maB = (Ind.bound   & (Pts.y1_kv >=  min(y1Plot)) & (Pts.y1_kv <=  max(y1Plot)) & (Pts.y2_kv <=  max(y2Plot)) & (Pts.y2_kv >=  min(y2Plot)));
    maC = (Ind.corners & (Pts.y1_kv >=  min(y1Plot)) & (Pts.y1_kv <=  max(y1Plot)) & (Pts.y2_kv <=  max(y2Plot)) & (Pts.y2_kv >=  min(y2Plot)));       
    
    rho_ic = rho_t(:,1);

    fl_norm = 0.1*max(max(max(abs(flux_t)))); xs = min(xl); ys = min(yl);
    
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
    
    y1_s  = Pts.y1_kv;  y2_s  = Pts.y2_kv;               
    
    pl_j = 1;
    for j = 1:d:n
                
        t      = plotTimes(j);        
        
        subplot(rows,cols,pl_j);        
        
        fl    = flux_t(:,1,j);
        fl_y1 = fl(1:N1*N2);
        fl_y2 = fl(N1*N2+1:end);                
        NormQuiverPlot(y1_s(ma3),y2_s(ma3),fl_y1(ma3),fl_y2(ma3),fl_norm,xs,ys,'b'); hold on;
        NormQuiverPlot(y1_s(maB),y2_s(maB),fl_y1(maB),fl_y2(maB),fl_norm,xs,ys,2,'k'); hold on;
        
        fl    = flux_t(:,2,j);
        fl_y1 = fl(1:N1*N2);
        fl_y2 = fl(N1*N2+1:end);                
        NormQuiverPlot(y1_s(ma3),y2_s(ma3),fl_y1(ma3),fl_y2(ma3),fl_norm,xs,ys,'r'); hold on;
        NormQuiverPlot(y1_s(maB),y2_s(maB),fl_y1(maB),fl_y2(maB),fl_norm,xs,ys,2,'k'); hold on;
        
                    
    
        y1M     = reshape(Interp.pts1,Interp.Nplot2,Interp.Nplot1);
        y2M     = reshape(Interp.pts2,Interp.Nplot2,Interp.Nplot1);            
    
        rho   = rho_t(:,1,j); z       = real(Interp.InterPol*rho);        
        [C,h]   = contour(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1),[0.1 0.2 0.3 0.4 0.5 0.6],'LineWidth',1.0,'color','b');    hold on;
        
        rho   = rho_t(:,2,j); z       = real(Interp.InterPol*rho);        
        [C,h]   = contour(y1M,y2M,reshape(z,Interp.Nplot2,Interp.Nplot1),[0.1 0.2 0.3 0.4 0.5 0.6],'LineWidth',1.0,'color','r');    hold on;
        

        %clabel(C,h,'FontSize',14);    

        
        pbaspect([(xl(2)-xl(1)) (yl(2)-yl(1)) 1]);
        xlim(xl); ylim(yl);
        xlabel('y_1');        ylabel('y_2');
        title(['t = ',num2str(t)]);
        pl_j = pl_j +1;
        
    end
        
    
end