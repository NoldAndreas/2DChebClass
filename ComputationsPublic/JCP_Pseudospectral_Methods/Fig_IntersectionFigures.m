function Fig_IntersectionFigures

    FigPlot(0,0.5,'InfAnnulus'); xlim([-2.5 2.5]); ylim([0 4]); SaveFigure('Fig_Intersection_1');
    FigPlot(0,1.1,'InfAnnulus'); xlim([-2.5 2.5]); ylim([0 4]); SaveFigure('Fig_Intersection_2');

    FigPlot(0,0.5,'Annulus'); xlim([-2.5 2.5]); ylim([0 4]); SaveFigure('Fig_Intersection_3');
    FigPlot(0,1.1,'Annulus'); xlim([-2.5 2.5]); ylim([0 4]); SaveFigure('Fig_Intersection_4');

    function FigPlot(y10,y20,shapeType)
       
        close all;

        L      = 1;
        R      = 1;
        bottom = 0;
        

        figure('color','white','Position',[0 0 800 500]);
        shape = struct('y2Min',bottom,'N',[10,10],'L1',1,'L2',2,'alpha',pi/2);
        HS    = HalfSpaceSkewed(shape);
        %HS    = HalfSpace(shape);
     %   HS.PlotGridLines();
       % plot([-10,10],[bottom,bottom],'b--','LineWidth',2);        

        Origin       = [y10,y20];
        N            = [15,16];
        sphere       = true;        
        
        shapeDC      = v2struct(Origin,N,sphere,L);
        shapeDC.RMin = 1;
        shapeDC.RMax = 2.5;
        shape        = str2func(shapeType);
        DC           = shape(shapeDC);   

        hold on
        area = Intersect(HS,DC);

        disp(['Area is: ',num2str(area.area),' with error from integration vector ',num2str(sum(area.int)-area.area)]);
        %disp(['Area from Int is: ',num2str(sum(area.int))]);
        
        for i = 1:length(area.shape.SubShape)
           area.shape.SubShape{i}.PlotGridLines(); 
           area.shape.SubShape{i}.PlotGrid(); 
        end

%        if(~isempty(area))                       
%            scatter(area.pts.y1_kv,area.pts.y2_kv-bottom,'ok');
%            %DC.PlotGrid();
%        end
        
        xlabel('$y_1$','Interpreter','Latex','fontsize',20);
        ylabel('$y_2$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',20);    
        pbaspect([1 1 1]);
        box on;

        axis equal
    end


end