clear all
close all

left  = -1;
right = 10;
top = 9;
bottom = 2;

figure
hold on

R = 2;

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


% SW
% y10 = left - 0.32*R;
% y20 = bottom - 0.5*R;

% NW
% y10 = left - 0.2*R;
% y20 = top + 0.5*R;

% NE
% y10 = right + 0.17*R;
% y20 = top + 0.4*R;

% SE
% y10 = right + 0.17*R;
% y20 = bottom - 0.1*R;

% SWin
% y10 = left + 0.32*R;
% y20 = bottom + 0.5*R;

% NWin
% y10 = left + 0.2*R;
% y20 = top - 0.5*R;

% NEin
% y10 = right - 0.17*R;
% y20 = top - 0.4*R;

% SEin
% y10 = right - 0.17*R;
% y20 = bottom + 0.1*R;

% SW double cut
% y10 = left + 0.95*R;
% y20 = bottom + 0.92*R;

% NW double cut
% y10 = left + 0.95*R;
% y20 = top - 0.92*R;

% NE double cut
% y10 = right - 0.95*R;
% y20 = top - 0.92*R;

% SE double cut
% y10 = right - 0.95*R;
% y20 = bottom + 0.92*R;

while(true)

[y10,y20] = ginput(1);

scatter(y10,y20,'b');
hold on
centralX = ( (y10 >= left + R) && (y10 <= right-R) );
centralY = ( (y20 >= bottom + R) && (y20 <= top-R) );
topY     = (y20 <= top + R) && (y20 >= top - R);
bottomY  = (y20 <= bottom + R) && (y20 >= bottom - R);
leftX    = (y10 <= left + R) && (y10 >= left - R);
rightX  = (y10 <= right + R) && (y10 >= right - R);

topIn  = (y20 <= top) && (y20 > top - R);
topOut = (y20 <= top + R) && (y20 > top); 

bottomIn  = (y20 < bottom + R) && (y20 >= bottom);
bottomOut = (y20 < bottom ) && (y20 >= bottom - R);

leftIn   = (y10 < left + R) && (y10 >= left);
leftOut  = (y10 < left) && (y10 >= left - R);

rightIn  = (y10 <= right) && (y10 > right - R);
rightOut = (y10 <= right + R) && (y10 > right);

containsNW = ( (y10-left)^2 + (y20-top)^2 < R^2);
containsNE = ( (y10-right)^2 + (y20-top)^2 < R^2);
containsSE = ( (y10-right)^2 + (y20-bottom)^2 < R^2);
containsSW = ( (y10-left)^2 + (y20-bottom)^2 < R^2);


Nin  = topIn && centralX;
Nout = topOut && (centralX || (leftIn && ~containsNW) || (rightIn && ~containsNE) );

Sin  = bottomIn && centralX;
Sout = bottomOut && (centralX || (leftIn &&~containsSW) || (rightIn && ~containsSE) );

Ein  = rightIn && centralY;
Eout = rightOut && (centralY || (topIn && ~containsNE) || (bottomIn && ~containsSE) );

Win  = leftIn && centralY;
Wout = leftOut && (centralY || (topIn && ~containsNW) || (bottomIn && ~containsSW) );


central = centralX && centralY;

NWw = containsNW;
NEw = containsNE;
SWw = containsSW;
SEw = containsSE;

NWwo = leftIn && topIn && ~containsNW;
NEwo = rightIn && topIn && ~containsNE;
SWwo = leftIn && bottomIn && ~containsSW;
SEwo = rightIn && bottomIn && ~containsSE;

corner = (NWw || NEw || SEw || SWw);

doubleCut = (NWwo || NEwo || SEwo || SWwo);

edgeIn = (Nin || Sin || Ein || Win);
edgeOut = (Nout || Sout || Eout || Wout);

if(central)
    shape.N = 20;
    shape.R = R;
    line = Circle(shape);
    line.Pts       = Pol2CartPts(line.Pts); 
    line.Pts.y1_kv = line.Pts.y1_kv + y10;
    line.Pts.y2_kv = line.Pts.y2_kv + y20;
    line.polar = 'cart';
elseif(corner)
    shape.N = 20;
    shape.Origin = [y10;y20];
    shape.R = R;
    
    if(SWw)
        shape.Corner = 'SW';
        shape.CornerPos = [left;bottom];
    elseif(NWw)
        shape.Corner = 'NW';
        shape.CornerPos = [left;top];
    elseif(NEw)
        shape.Corner = 'NE';
        shape.CornerPos = [right;top];
    elseif(SEw)
        shape.Corner = 'SE';
        shape.CornerPos = [right;bottom];
    end
    
    line = CornerCircle(shape);
    
elseif(doubleCut)
    shape.N = 20;
    shape.Origin = [y10;y20];
    shape.R = R;
    
    if(SWwo)
        shape.Corner = 'SW';
        shape.CornerPos = [left;bottom];
    elseif(NWwo)
        shape.Corner = 'NW';
        shape.CornerPos = [left;top];
    elseif(NEwo)
        shape.Corner = 'NE';
        shape.CornerPos = [right;top];
    elseif(SEwo)
        shape.Corner = 'SE';
        shape.CornerPos = [right;bottom];
    end
    
    line = DoubleCutCircle(shape);

elseif(edgeIn || edgeOut)
    shape.N = 20;
    shape.R  = R;
    shape.Origin = [y10;y20];
    if(Sin || Sout)
        shape.h = y20 - bottom;
        shape.WallPos = 'S';
    elseif(Nin || Nout)
        shape.h = top - y20;
        shape.WallPos = 'N';
    elseif(Ein || Eout)
        shape.h = right - y10;
        shape.WallPos = 'E';
    elseif(Win || Wout)
        shape.h = y10 - left;
        shape.WallPos = 'W';
    end
    
    line = ArcNSEW(shape); 

else
    
    line = [];
end

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