clear all; close all;

R      = 2;
bottom = 2;

shape = struct('y2Min',bottom,'N',[10,10],'L1',1,'L2',2);
HS    = HalfSpace(shape);
plot([-6,6],[bottom,bottom],'r--','LineWidth',2);

while(true)
    xlim([-3,3])
    ylim([bottom-R,bottom+3*R])
    axis equal

    
    [y10,y20] = ginput(1);

    Origin = [y10,y20];
    N      = [10,10];
    DC     = Disc(v2struct(Origin,R,N));   

    hold on
    area = Intersect(HS,DC,[]);
    
    disp(['Area is: ',num2str(area.area)]);
    disp(['Area from Int is: ',num2str(sum(area.int))]);

    if(~isempty(area))                       
        scatter(area.pts.y1_kv,area.pts.y2_kv,'r');
        %DC.PlotGrid();
    end

end