clear all
close all

R = 2;

xw = 1;
yw = 3;

x0 = xw + 0.9*R;
y0 = yw + 0.8*R;

xip = x0 + sqrt(R^2 - (y0-yw)^2);
yip = y0 + sqrt(R^2 - (x0-xw)^2);

xim = x0 - sqrt(R^2 - (y0-yw)^2);
yim = y0 - sqrt(R^2 - (x0-xw)^2);


figure
hold on


geomD.N = [20;20];
geomD.Origin = [x0;y0];
geomD.R = R;

D = Disc(geomD);
D.Pts = Pol2CartPts(D.Pts);
D.polar = 'cart';
D.Pts.y1_kv = D.Pts.y1_kv + x0;
D.Pts.y2_kv = D.Pts.y2_kv + y0;

%D.PlotGrid;

scatter([xip, xim],[yw,yw],'MarkerFaceColor','r');
scatter([xw,xw],[yip,yim],'MarkerFaceColor','b');


plot([xw xw],ylim(gca),'r')
plot(xlim(gca),[yw,yw],'r')

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

