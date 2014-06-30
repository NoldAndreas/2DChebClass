clear all; close all;

R      = 2;
left   = -1; right = 10;
bottom = 2;  top   = 9;

shape = struct('y1Min',left,'y1Max',right,...
               'y2Min',bottom,'y2Max',top,...
               'N',[10,10]);

BX = Box(shape);
BX.SetUpBorders(20);

figure;
BX.PlotBorders();

axis equal
xlim([left-R,right+R])
ylim([bottom-R,top+R])


plot([left+R left+R],ylim(gca),'Color','r','LineWidth',2,'LineStyle','--')
plot([right-R right-R],ylim(gca),'Color','r','LineWidth',2,'LineStyle','--')
plot(xlim(gca),[top-R,top-R],'Color','r','LineWidth',2,'LineStyle','--')
plot(xlim(gca),[bottom+R,bottom+R],'Color','r','LineWidth',2,'LineStyle','--')

while(true)

    [y10,y20] = ginput(1);

    Origin = [y10,y20];
    N      = [10,10];
    DC     = Disc(v2struct(Origin,R,N));   

    hold on
    area = Intersect(BX,DC,[]);
    
    disp(['Area is: ',num2str(area.area)]);

    if(~isempty(area))                       
        scatter(area.pts.y1_kv,area.pts.y2_kv,'r');
        %DC.PlotGrid();
    end

end