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

scatter(y10,y20,'b');
hold on



shapeC.N = 20;
shapeC.R = R;
lineC = Circle(shapeC);

if(~isempty(line))
    
    dataCircle.pts = line.Pts;
    
    fullCircle.pts = Pol2CartPts(lineC.Pts);
    fullCircle.pts.y1_kv = fullCircle.pts.y1_kv + y10;
    fullCircle.pts.y2_kv = fullCircle.pts.y2_kv + y20;
    
    scatter(fullCircle.pts.y1_kv,fullCircle.pts.y2_kv,'b');
    hold on;
    scatter(dataCircle.pts.y1_kv,dataCircle.pts.y2_kv,'r');
    
end

end