clear all
close all
figure

Phys_Area = struct('shape','HalfSpace_FMT','N',[10;10],'L1',2,'L2',2, ...
                   'y2wall',0,'N2bound',10,'h',1,'L2_AD',1);

Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                   'y2Min',0,'y2Max',5,'N2',100);

R = 0.5;
               
IDC = HalfSpace_FMT(Phys_Area,R);
    
[Pts,Diff,Int,Ind,~] = IDC.ComputeAll(Plot_Area);

IDC.PlotGrid;
xlim([-15,15]); ylim([-10,20]);

%[x0,y0] = ginput(1);

x0 = 0;
y0 = 0.8848;

opts.offset_y2 = y0;

%% Annulus - 1
annGeom.N      = [10;10];
annGeom.RMin   = R;
annGeom.L      = 2;
%annGeom.Origin = [x0;y0];

fullArea       = InfAnnulus(annGeom);
cutArea        = Intersect(IDC,fullArea,opts);

hold on
annPts         = fullArea.GetCartPts();
scatter(annPts.y1_kv+x0,annPts.y2_kv+y0,'b')

scatter(cutArea.pts.y1_kv+x0,cutArea.pts.y2_kv,'r')

% [cutArea.pts.y1_kv cutArea.pts.y2_kv]
% 
% pause

figure;
IDC.PlotGrid;
xlim([-15,15]); ylim([-10,20]);

%% Annulus - 2
annGeom.Origin = [x0;y0];

fullArea       = InfAnnulus(annGeom);
cutArea        = Intersect(IDC,fullArea,[]);

hold on
annPts         = fullArea.GetCartPts();
scatter(annPts.y1_kv,annPts.y2_kv,'b')

scatter(cutArea.pts.y1_kv,cutArea.pts.y2_kv,'r')

%% Disc
discGeom.N = [20;20];
discGeom.R = R/2;
%annGeom.Origin = [x0;y0];

discArea    = Disc(discGeom);
cutDiscArea = Intersect(IDC,discArea,opts);

discPts     = discArea.GetCartPts();
scatter(discPts.y1_kv + x0,discPts.y2_kv+y0,'g')

scatter(cutDiscArea.pts.y1_kv + x0,cutDiscArea.pts.y2_kv,'m')