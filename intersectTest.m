clear all
close all

R = 2;

xw = 1;
yw = 3;

x0 = xw + R*(2*rand-1);
y0 = yw + R*(2*rand-1);

xi = x0 + sqrt(R^2 - (y0-yw)^2);
yi = y0 + sqrt(R^2 - (x0-xw)^2);

xp = 1/2*(xi+xw);
yp = 1/2*(yi+yw);

d0p = sqrt( (x0- xp)^2 + (y0 - yp)^2 );

NT = [10;10];

YT = [xw yw; xw yi; xi yw];

geomT.Y = YT;

t = Triangle(geomT,NT);

figure
t.PlotGrid;
hold all

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

geomD.N = [20;20];
geomD.Origin = [x0;y0];
geomD.R = R;

D = Disc(geomD);
D.Pts = Pol2CartPts(D.Pts);
D.polar = 'cart';
D.Pts.y1_kv = D.Pts.y1_kv + x0;
D.Pts.y2_kv = D.Pts.y2_kv + y0;

D.PlotGrid;



plot([xw xw],ylim(gca),'r')
plot(xlim(gca),[yw,yw],'r')