clear all; close all;

R_in   = 1;
R_out  = 3;
R      = 3;%0.5;
bottom = 0;

figure('color','white','Position',[0 0 800 500]);
shape = struct('y2Min',bottom,'N',[10,10],'L1',1,'L2',2,'alpha',pi/4);
HS    = HalfSpaceSkewed(shape);
%HS    = HalfSpace(shape);
HS.PlotGridLines();

plot([-10,10],[bottom,bottom],'b--','LineWidth',2);


xlim([-5,5])
ylim([-1.3 2])
ylim([bottom-R-1,bottom+3*R])
xlabel('$x$','Interpreter','Latex','fontsize',20);
ylabel('$y$','Interpreter','Latex','fontsize',20);

%axis equal
pbaspect([10 4 1]);
set(gca,'fontsize',20);


[y10,y20] = ginput(1);   
while(y20 > bottom - R)
    
    Origin = [y10,y20];
    N      = [10,10];
    sphere = false;
    
    DC     = Annulus(v2struct(Origin,R_in,R_out,N));       
    %DC     = Disc(v2struct(Origin,R,N,sphere));       
    %DC     = Sphere(v2struct(Origin,R,N,sphere));       
    %theta1 = 0;  theta2 = pi;    
    %DC           = Ball(v2struct(Origin,N,sphere,theta1,theta2,R));       

    hold on
    area = Intersect(HS,DC);
    
    disp(['Area is: ',num2str(area.area),' with error from integration vector ',num2str(sum(area.int)-area.area)]);
    %disp(['Area from Int is: ',num2str(sum(area.int))]);

    if(~isempty(area))                       
        scatter(area.pts.y1_kv,area.pts.y2_kv,'r');
        %DC.PlotGrid();
    end
    
	xlim([-3,3])
    ylim([-1.5 1.5])
    ylim([bottom-R,bottom+3*R])
    xlabel('$x/\sigma$','Interpreter','Latex','fontsize',20);
    ylabel('$y/\sigma$','Interpreter','Latex','fontsize',20);
    set(gca,'fontsize',20);
    pbaspect([6 4 1]);
    %axis equal

    [y10,y20] = ginput(1);   
end
ylim([-1 3]);
print2eps(['HalfSpaceVsDisc_CollocationPoints'],gcf);        
saveas(gcf,['HalfSpaceVsDisc_CollocationPoints.fig']);   