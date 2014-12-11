
close all;
global dirData
AddPaths();        
ChangeDirData([dirData filesep 'CahnHilliard_InnerRegion'],'ORG');    

PhysArea = struct('N',[50,40],'y2Min',0,'y2Max',20,...
                  'L1',7,'IntInterval',[-10,10]);%,'NBorder',[30,200,30,200]);

PlotArea = struct('y1Min',-15,'y1Max',15,'N1',80,'N2',80);   
SubArea  = struct('shape','Box','N',[60,60],...
                  'y1Min',-5,'y1Max',5,'y2Min',0,'y2Max',10);   	

optsNum  = v2struct(PhysArea,PlotArea);                   	

optsPhys = struct('thetaEq',pi/2,...       
                   'theta',90*pi/180,...
                   'Cak',0.01,'Cn',1,...
                   'UWall',1,...                       
                   'mobility',10,...
                   'nParticles',0);

config = v2struct(optsPhys,optsNum);   

opts = struct('noIterations',20,'lambda',0.8,'Seppecher_red',1);


for y2Max = 5%20%10:5:15
    
    config.optsNum.PhysArea.y2Max = y2Max;
    
    DI = DiffuseInterfaceBinaryFluid(config);
    DI.Preprocess();
    
    DI.IterationStepFullProblem(opts);    

    opts.Seppecher_red = 2;
    opts.lambda        = 0.6;    
    %DI.optsPhys.mobility = m;
    DI.IterationStepFullProblem(opts);    
    DI.PlotResults();	  
    
    clear('DI');
end
% 
% opts.solveSquare   = false;
% DI.IterationStepFullProblem(opts);    

DI.optsNum.SubArea = SubArea;
DI.Preprocess_SubArea();
DI.PostProcess_Flux;
DI.IDC.SetUpBorders([30,1000,30,200]);
DI.FindAB();

%DI.FindStagnationPoint();
DI.PlotResults();	               
