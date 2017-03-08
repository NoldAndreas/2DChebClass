function StreamlinePlotInterp_NSpecies(Interp,Pts,flux,lw,c)
    
    nSpecies=size(flux,2);
    
    y1_s= Interp.pts1; %Pts.y1_kv;
    y2_s= Interp.pts2; %Pts.y2_kv;
    
    N1 = Pts.N1;    N2 = Pts.N2;
    
    fl_y1          = Interp.InterPol*flux(1:N1*N2,:);
    fl_y2          = Interp.InterPol*flux(N1*N2+1:end,:);
    
    for iSpecies=1:nSpecies
         
        fl_y1S=fl_y1(:,iSpecies);
        fl_y2S=fl_y2(:,iSpecies);
        
        if(nargin == 7) 
            quiver([pos_norm(1);y1_s],[pos_norm(2);y2_s],...
                [fl_norm(1);fl_y1S],[fl_norm(2);fl_y2S],...
                'LineWidth',lw,'Color',c{iSpecies});
        else
            quiver([pos_norm(1);y1_s],[pos_norm(2);y2_s],...
                [fl_norm(1);fl_y1S],[fl_norm(2);fl_y2S]);
        end
        
        hold on;
        
    end
    
    xl = [min(y1_s) max(y1_s)];
    yl = [min(y2_s) max(y2_s)];
    
    xlim(xl);
    ylim(yl);
    
    pbaspect([(xl(2)-xl(1))  (yl(2)-yl(1)) 1]);
    xlabel('y_1');
    ylabel('y_2');
end