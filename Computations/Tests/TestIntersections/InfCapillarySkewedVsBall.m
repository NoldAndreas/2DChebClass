clear all; close all;

R      = 1;
bottom = 1;
top    = 5;

IntersectShape = 'Ball'; %Disc

figure('color','white','Position',[0 0 800 500]);
shape = struct('y2Min',bottom,'y2Max',top,'N',[6,6],'L1',1,'L2',2,...
               'alpha',70/180*pi);
           
HS    = InfCapillarySkewed(shape);
HS.PlotGridLines();
plot([-10,10],[bottom,bottom],'b--','LineWidth',2); hold on;
plot([-10,10],[top,top],'b--','LineWidth',2);

xlim([-5,5])
ylim([(bottom-R-1) (top + R + 1)])
xlabel('$x$','Interpreter','Latex','fontsize',20);
ylabel('$y$','Interpreter','Latex','fontsize',20);
set(gca,'fontsize',20);
pbaspect([10 (top-bottom+2*R+2) 1]);

[y10,y20] = ginput(1);   
while(y20 > bottom - R)
    
    Origin       = [y10,y20];
    N            = [10,10];
    sphere       = true;
    
    if(strcmp(IntersectShape,'Disc'))
        DC     = Disc(v2struct(Origin,R,N,sphere));       
    elseif(strcmp(IntersectShape,'Ball'))    
        theta1 = 0;  theta2 = pi;    
        DC           = Ball(v2struct(Origin,N,sphere,theta1,theta2,R));   
	elseif(strcmp(IntersectShape,'InfAnnulus')) 
        L            = 1;
        shapeDC      = v2struct(Origin,N,sphere,L);
        shapeDC.RMin = R;
        DC           = InfAnnulus(shapeDC);  
    end
        
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
    pbaspect([10 (top-bottom+2*R+2) 1]);

    [y10,y20] = ginput(1);   
end
xlim([-5 5]);
ylim([-1 9]);
axis equal
print2eps(['HalfSpaceVsInfAnnulus_CollocationPoints'],gcf);        
saveas(gcf,['HalfSpaceVsInfAnnulus_CollocationPoints.fig']);   