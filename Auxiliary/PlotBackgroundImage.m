function PlotBackgroundImage(figFile,shift)
     background=imread(figFile);
          
     %temp = shift.temp/2;
     xmin = shift.xmin; %temp
     xmax = shift.xmax;
     ymin = shift.ymin;
     ymax = shift.ymax;

     % flip image

     imagesc([xmin xmax], [ymin ymax], flipud(background));
     colormap(gray);

     set(gca,'ydir','normal');

     hold on
     
	 plot([shift.xref shift.xref],[ymin,ymax],'g','linewidth',1.5);
     plot([xmin xmax],[shift.yref shift.yref],'g','linewidth',1.5);
    
end