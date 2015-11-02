close all;

L      = 1;
R      = 1;
bottom = 0.5;

figure('color','white','Position',[0 0 800 500]);
shape = struct('y2Min',bottom,'N',[10,10],'L1',1,'L2',2,'alpha',pi/4);
HS    = HalfSpaceSkewed(shape);
%HS    = HalfSpace(shape);
HS.PlotGridLines();
plot([-10,10],[bottom,bottom],'b--','LineWidth',2);

xlim([-5,5])
ylim([-1 9])
xlabel('$x$','Interpreter','Latex','fontsize',20);
ylabel('$y$','Interpreter','Latex','fontsize',20);
set(gca,'fontsize',20);
pbaspect([1 1 1]);

[y10,y20] = ginput(1);   
while(y20 > bottom - R)
    
    Origin       = [y10,y20];
    N            = [15,16];
    sphere       = true;
    shapeDC      = v2struct(Origin,N,sphere,L);
    shapeDC.RMin = R;
    DC           = InfAnnulus(shapeDC);   

    hold on
    area = Intersect(HS,DC);
    
    disp(['Area is: ',num2str(area.area),' with error from integration vector ',num2str(sum(area.int)-area.area)]);
    %disp(['Area from Int is: ',num2str(sum(area.int))]);

    if(~isempty(area))                       
        scatter(area.pts.y1_kv,area.pts.y2_kv,'r');
        %DC.PlotGrid();
    end
    
	xlim([-5,5])
    ylim([-1 9])    
    xlabel('$x$','Interpreter','Latex','fontsize',20);
    ylabel('$y$','Interpreter','Latex','fontsize',20);
    set(gca,'fontsize',20);    
    pbaspect([1 1 1]);

    [y10,y20] = ginput(1);   
end
xlim([-5 5]);
ylim([-1 9]);
axis equal
print2eps(['HalfSpaceVsInfAnnulus_CollocationPoints'],gcf);        
saveas(gcf,['HalfSpaceVsInfAnnulus_CollocationPoints.fig']);   