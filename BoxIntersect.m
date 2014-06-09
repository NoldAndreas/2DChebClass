clear all
close all

left  = -1;
right = 10;
top = 8;
bottom = 2;

figure
hold on

R = 1.5;

axis equal
xlim([left-R,right+R])
ylim([bottom-R,top+R])

plot([left right],[top top],'Color','r','LineWidth',2)
plot([left right],[bottom bottom],'Color','r','LineWidth',2)
plot([left left],[top bottom],'Color','r','LineWidth',2)
plot([right right],[top bottom],'Color','r','LineWidth',2)

plot([left+R left+R],ylim(gca),'Color','b','LineWidth',2,'LineStyle','--')
plot([right-R right-R],ylim(gca),'Color','b','LineWidth',2,'LineStyle','--')
plot(xlim(gca),[top-R,top-R],'Color','b','LineWidth',2,'LineStyle','--')
plot(xlim(gca),[bottom+R,bottom+R],'Color','b','LineWidth',2,'LineStyle','--')



while(true)
    [x0,y0] = ginput(1);
    scatter(x0,y0,'r');
    centralX = ( (x0 >= left + R) && (x0 <= right-R) );
    centralY = ( (y0 >= bottom + R) && (y0 <= top-R) );
    topY     = (y0 <= top + R) && (y0 >= top - R);
    bottomY  = (y0 <= bottom + R) && (y0 >= bottom - R);
    leftX    = (x0 <= left + R) && (x0 >= left - R);
    rightX  = (x0 <= right + R) && (x0 >= right - R);

    central = centralX && centralY;
    N = centralX && topY;
    S = centralX && bottomY;
    E = rightX & centralY;
    W = leftX && centralY;
    
    NW = leftX && topY && ( (x0-left)^2 + (y0-top)^2 < R^2);
    NE = rightX && topY && ( (x0-right)^2 + (y0-top)^2 < R^2);
    SE = rightX && bottomY && ( (x0-right)^2 + (y0-bottom)^2 < R^2);
    SW = leftX && bottomY && ( (x0-left)^2 + (y0-bottom)^2 < R^2);
    
end