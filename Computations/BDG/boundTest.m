PhysArea = struct('y1Min',-inf,'y1Max',inf,'L1',3,'N',[12;10],...
                       'y2Min',-2,'y2Max',2);
                   
PlotArea = struct('y1Min',-5,'y1Max',5,'y2Min',-2,'y2Max',2, ...
                       'N1',100,'N2',100);
    
FexNum = struct('Fex','Meanfield','N',[20,20],'L1',6,'L2',1);

shape        = PhysArea;
shape.Conv   = FexNum;

IDC                   = InfCapillary(shape);

[Pts,Diff,Int,Ind,~]  = IDC.ComputeAll(PlotArea);

%Ind.finite  = Ind.bound & (isfinite(Pts.y1_kv) | isfinite(Pts.y2_kv));
% Ind.infinite  = Ind.bound & (~isfinite(Pts.y1_kv) | ~isfinite(Pts.y2_kv));
% 
% Ind.finite   = Ind.bound & ~Ind.infinite;
% 
% [Ind.finite (Ind.top | Ind.bottom) Ind.infinite (Ind.left | Ind.right)];

% Ind.normal
% Ind.normalTop
% Ind.normalBottom
% Ind.normalLeft
% Ind.normalRight

temp = (1:2*PhysArea.N(1)*PhysArea.N(2))';

% Ind.normalTop*temp 
% Ind.normalBottom*temp 
% Ind.normalLeft*temp 
% Ind.normalRight*temp


% size(Ind.bound)
% size(Ind.finite)
% nnz(Ind.finite)

[Ind.finite (Ind.top | Ind.bottom) Ind.infinite (Ind.left | Ind.right)];