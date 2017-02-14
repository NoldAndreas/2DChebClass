function JCP_639_2017_Fig4_IntersectionFigures

    AddPaths('JCP_639_2017');
    close all;

    FigPlot(0,0.5,'InfAnnulus'); xlim([-2.5 2.5]); ylim([0 4]); SaveFigure('JCP_639_2017_Fig4_IntersectionFigures_a');
    FigPlot(0,1.1,'InfAnnulus'); xlim([-2.5 2.5]); ylim([0 4]); SaveFigure('JCP_639_2017_Fig4_IntersectionFigures_b');

    FigPlot(0,0.5,'Annulus'); xlim([-2.5 2.5]); ylim([0 4]); SaveFigure('JCP_639_2017_Fig4_IntersectionFigures_c');
    FigPlot(0,1.1,'Annulus'); xlim([-2.5 2.5]); ylim([0 4]); SaveFigure('JCP_639_2017_Fig4_IntersectionFigures_d');

    function FigPlot(y10,y20,shapeType)
       
        shape        = struct('y2Min',0,'N',[10,10],'L1',1,'L2',2,'alpha',pi/2);
        HS           = HalfSpaceSkewed(shape);
        
        shapeDC      = struct('N',[15,16],...
                              'Origin',[y10,y20],...
                              'sphere',true,...
                              'L',1,...
                              'RMin',1,'RMax',2.5);        
        objType      = str2func(shapeType);
        DC           = objType(shapeDC);   

        area         = Intersect(HS,DC);

        %**** Plotting ****
        figure('color','white','Position',[0 0 800 500]);
        for i = 1:length(area.shape.SubShape)
           area.shape.SubShape{i}.PlotGridLines(); 
           area.shape.SubShape{i}.PlotGrid(); 
        end        
        xlabel('$y_1$','Interpreter','Latex','fontsize',20);
        ylabel('$y_2$','Interpreter','Latex','fontsize',20);
        set(gca,'fontsize',20);    
        pbaspect([1 1 1]);
        box on; axis equal;
    end

end