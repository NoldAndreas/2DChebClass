clear all
close all
figure

Phys_Area = struct('shape','HalfSpace_FMT','N',[30;30],'L1',2,'L2',2, ...
                   'y2wall',0,'N2bound',10,'h',1,'L2_AD',1);

Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                   'y2Min',0.5,'y2Max',5,'N2',100);

R = 5;
               
IDC = HalfSpace_FMT(Phys_Area,R,1e-6);
    
[Pts,Diff,Int,Ind,~] = IDC.ComputeAll(Plot_Area);

IDC.PlotGrid;

xlim([-15,15]);
ylim([-10,20]);

[x0,y0] = ginput(1);

annGeom.N = [20;20];
annGeom.RMin = R;
annGeom.L = 2;
annGeom.Origin = [x0;y0-R];

fullArea = InfAnnulus(annGeom);

opts.offset_y2 = R;

cutArea = Intersect(IDC,fullArea,opts);

hold on

annPts = Pol2CartPts(fullArea.Pts);
annPts.y1_kv = annPts.y1_kv + x0;
annPts.y2_kv = annPts.y2_kv + y0;
scatter(annPts.y1_kv,annPts.y2_kv,'b')

scatter(cutArea.pts.y1_kv,cutArea.pts.y2_kv,'r')


discGeom.N = [20;20];
discGeom.R = R/2;
annGeom.Origin = [x0;y0];

discArea = Disc(discGeom);

opts.offset_y2 = y0;

cutDiscArea = Intersect(IDC,discArea,opts);


discPts = Pol2CartPts(discArea.Pts);
discPts.y1_kv = discPts.y1_kv + x0;
discPts.y2_kv = discPts.y2_kv + y0;
scatter(discPts.y1_kv,discPts.y2_kv,'g')

scatter(cutDiscArea.pts.y1_kv + x0,cutDiscArea.pts.y2_kv,'m')
