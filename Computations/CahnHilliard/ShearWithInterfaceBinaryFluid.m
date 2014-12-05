    close all;
    
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_ShearedInterface' filesep 'BinaryFluid'],'ORG');
    %% Parameters    
    PhysArea = struct('N',[50,40],'y2Min',0,'y2Max',[],'L1',10,'IntInterval',[-10 10]); %12,80,50                      
    PlotArea = struct('y1Min',-20,'y1Max',20,'N1',80,'N2',80);   
    SubArea  = struct('shape','Box','N',[60,60],...
                      'y1Min',-1,'y1Max',1,...
                      'y2Min',4,'y2Max',6);   
	optsNum  = v2struct(PhysArea,PlotArea);
    
    optsPhys = struct('thetaEq',[pi/2,pi/2],...
                      'eta',1,...
                      'Cak',0.03,'Cn',1,...
                      'UWall',[1,-1],...                      .
                      'mobility',10,...
                      'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
                      
    for y2M = 20%:2:20%22:2:30%15:2.5:30
        
        config.optsNum.PhysArea.y2Max = y2M;        
        
        DI = DiffuseInterfaceBinaryFluid(config);
        DI.Preprocess();               
        
        opts = struct('noIterations',20,'lambda',0.8,'solveSquare',true);
        DI.IterationStepFullProblem(opts);    
	
        opts = struct('noIterations',20,'lambda',0.8,'solveSquare',false);
        DI.IterationStepFullProblem(opts);    
                
        DI.optsNum.SubArea = SubArea;
        DI.Preprocess_SubArea();

        DI.PlotResults();    
        clear 'DI'
    end



