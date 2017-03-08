function PlotInterfaceAnalysisY1(this)

    filmThickness = this.filmThickness;
    isoInterface  = this.IsolineInterface;
    y1            = this.y1; 
    
    %********************************************************************
    %Film thickness            
    IntM = zeros(length(y1));
	vh   = zeros(1,length(y1));
        
    for i = 2:length(y1)        
        %mark_II  = (y1>3);
        %y1_II    = y1(mark_II);
        %Diff_II  = barychebdiff(y1_II,2);
        %D_II     = Diff_II.Dx; D2_II = Diff_II.DDx;        

        %Matrix for integral int(v(y),y=y1..Y)

        
        h           = y1(i) - y1(i-1);
        vh([i-1,i]) = vh([i-1,i]) + h/2;
        IntM(i,:)   = vh;
    end
    if(~isempty(filmThickness))
        figure('Color','white','Position',[0 0 600 1000],'name','Potentials');    
        ellP = (D*filmThickness)./((1+(D*filmThickness).^2).^0.5);
        ell  = IntM*ellP + filmThickness(1);


        subplot(1,2,1); 
        plot(y1,filmThickness); hold on;  plot(y1,ell,'m'); %plot(y1,filmThickness_II,'r'); 
        plot(y1,isoInterface,'r');

        subplot(1,2,2); 
        plot(y1,atan(this.DiffY1*filmThickness)*180/pi); hold on; % plot(y1(mark_II),D_II*filmThickness_II(mark_II),'r');
        plot(y1,atan(this.DiffY1*isoInterface)*180/pi,'r'); hold on; % plot(y1(mark_II),D_II*filmThickness_II(mark_II),'r');
        ylabel('Slope [deg]');

        %subplot(2,2,3); plot(y1,D2*filmThickness); hold on; % plot(y1(mark_II),D2_II*filmThickness_II(mark_II),'r');
        %plot(y1,D2*isoInterface);
        %subplot(2,2,4); plot(y1,(D*filmThickness)./((1+(D*filmThickness).^2).^0.5));    

        %subplot(2,2,4); plot(y1,ell);
    end
       
end
