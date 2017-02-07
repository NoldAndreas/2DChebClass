function StreamlinePlotInterp(Interp,Pts,flux,lw,c)
    

    y1_s= Interp.pts1; y2_s= Interp.pts2; 
    
    N1 = Pts.N1;    N2 = Pts.N2;
    
    fl_y1          = Interp.InterPol*flux(1:N1*N2,:);
    fl_y2          = Interp.InterPol*flux(N1*N2+1:end,:);
    
    y1M     = reshape(y1_s,Interp.Nplot2,Interp.Nplot1);
    y2M     = reshape(y2_s,Interp.Nplot2,Interp.Nplot1);
    fl_y1M  = reshape(fl_y1,Interp.Nplot2,Interp.Nplot1);
    fl_y2M  = reshape(fl_y2,Interp.Nplot2,Interp.Nplot1);
    
    %Define starting points for streamlines
    dy = 0.3;
    Ns = ceil((max(y2_s)-min(y2_s))/dy);
    Sx = ones(Ns,1)*max(y1_s);
    Sy = min(y2_s) + (max(y2_s)-min(y2_s))*(0:Ns-1)'/Ns;    
    
    Ns = ceil((max(y2_s)-min(y2_s))/dy);
    Sx = [Sx;ones(Ns,1)*min(y1_s)];
    Sy = [Sy;min(y2_s) + (max(y2_s)-min(y2_s))*(0:Ns-1)'/Ns];
    
    dx = 0.5;
    Ns = ceil((max(y1_s)-min(y1_s))/dx);
    Sx = [Sx;min(y1_s) + (max(y1_s)-min(y1_s))*(0:Ns-1)'/Ns];
    Sy = [Sy;max(y2_s)*ones(Ns,1)];
    
    
    
    if(nargin == 5)         
       streamline(y1M,y2M,fl_y1M,fl_y2M,Sx,Sy,'LineWidth',lw,'Color',c);
    else
       streamline(y1M,y2M,fl_y1M,fl_y2M,Sx,Sy);
   end

    
    %Configuring the plot
    xl = [min(y1_s) max(y1_s)];    yl = [min(y2_s) max(y2_s)];    
    xlim(xl);    ylim(yl);
    
    pbaspect([(xl(2)-xl(1))  (yl(2)-yl(1)) 1]);
    xlabel('y_1');  ylabel('y_2');
end