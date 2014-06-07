clear all
close all

R = 2;

xw = 1;
yw = 3;

x0 = xw + R*(2*rand-1);
y0 = yw + R*(2*rand-1);


figure
hold on

xlim([xw-2*R,xw+2*R])
ylim([yw-2*R,yw+2*R])

plot([xw xw],ylim(gca),'Color','r','LineWidth',2)
plot(xlim(gca),[yw,yw],'Color','r','LineWidth',2)

plot([xw-R xw-R],ylim(gca),'Color','b','LineWidth',2,'LineStyle','--')
plot([xw+R xw+R],ylim(gca),'Color','b','LineWidth',2,'LineStyle','--')
plot(xlim(gca),[yw-R,yw-R],'Color','b','LineWidth',2,'LineStyle','--')
plot(xlim(gca),[yw+R,yw+R],'Color','b','LineWidth',2,'LineStyle','--')

[x0,y0] = ginput(1);

originInDisc = ( (x0-xw)^2 + (y0-yw)^2 <= R^2);

if(originInDisc)

    xi = x0 + sqrt(R^2 - (y0-yw)^2);
    yi = y0 + sqrt(R^2 - (x0-xw)^2);

    xp = 1/2*(xi+xw);
    yp = 1/2*(yi+yw);

    d0p = sqrt( (x0- xp)^2 + (y0 - yp)^2 );

    NT = [10;10];

    YT = [xw yw; xw yi; xi yw];

    geomT.Y = YT;

    t = Triangle(geomT,NT);


    t.PlotGrid;

    h = d0p;

    geomS.N = [20;20];
    geomS.R = R;
    geomS.h = h;

    theta = -atan((yi-yw)/(xi-xw));

    s = Segment(geomS);

    y1_kv = s.Pts.y1_kv;
    y2_kv = s.Pts.y2_kv - h;

    y1_kvs = cos(theta)*y1_kv - sin(theta)*y2_kv + xp;
    y2_kvs = sin(theta)*y1_kv + cos(theta)*y2_kv + yp;

    s.Pts.y1_kv = y1_kvs;
    s.Pts.y2_kv = y2_kvs;

    s.PlotGrid;

else

    xip = x0 + sqrt(R^2 - (y0-yw)^2);
    yip = y0 + sqrt(R^2 - (x0-xw)^2);

    xim = x0 - sqrt(R^2 - (y0-yw)^2);
    yim = y0 - sqrt(R^2 - (x0-xw)^2);

    NT = [10;10];

    geomT1.Y = [xw yim; xw yip; x0 y0];

    t1 = Triangle(geomT1,NT);
    t1.PlotGrid;

    geomT2.Y = [xim yw; x0 y0; xip yw];

    t2 = Triangle(geomT2,NT);
    t2.PlotGrid;


    theta1 = pi + atan((y0-yw)/(x0-xim));
    theta2 = 3*pi/2 - atan((x0-xw)/(y0-yim));

    geomW1.N = [20;20];
    geomW1.R = R;
    geomW1.th1 = theta1;
    geomW1.th2 = theta2;

    w1 = Wedge(geomW1);

    w1.Pts  = Pol2CartPts(w1.Pts);
    w1.polar = 'cart';

    w1.Pts.y1_kv = w1.Pts.y1_kv + x0;
    w1.Pts.y2_kv = w1.Pts.y2_kv + y0;
    w1.PlotGrid;

    theta1 = -pi/2 + atan((xip-x0)/(y0-yw));
    theta2 = pi - atan((yip-y0)/(x0-xw));


    geomW2.N = [20;20];
    geomW2.R = R;
    geomW2.th1 = theta2;
    geomW2.th2 = theta1;

    w2 = Wedge(geomW2);

    w2.Pts  = Pol2CartPts(w2.Pts);
    w2.polar = 'cart';

    w2.Pts.y1_kv = w2.Pts.y1_kv + x0;
    w2.Pts.y2_kv = w2.Pts.y2_kv + y0;
    w2.PlotGrid;

end

geomD.N = [20;20];
geomD.Origin = [x0;y0];
geomD.R = R;

D = Disc(geomD);
D.Pts = Pol2CartPts(D.Pts);
D.polar = 'cart';
D.Pts.y1_kv = D.Pts.y1_kv + x0;
D.Pts.y2_kv = D.Pts.y2_kv + y0;

%D.PlotGrid;
   

