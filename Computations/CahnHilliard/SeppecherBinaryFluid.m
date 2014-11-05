function SeppecherBinaryFluid()

    close all;
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_InnerRegion'],'ORG');    
        
    PhysArea = struct('N',[50,40],'y2Min',0,'y2Max',10,...
                      'L1',7,'IntInterval',[-10,10]);%,'NBorder',[30,200,30,200]);

    PlotArea = struct('y1Min',-15,'y1Max',15,'N1',80,'N2',80);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('thetaEq',pi/2,...                         
                       'Cak',0.02,'Cn',1,...
                       'UWall',1,...                       
                       'mobility',10,...
                       'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
            
    DI = DiffuseInterfaceBinaryFluid(config);
    DI.Preprocess();
    
    opts = struct('noIterations',10,'lambda',0.5,'solveSquare',true);
    DI.IterationStepFullProblem(opts);    
    
    opts = struct('noIterations',5,'lambda',0.8,'solveSquare',true);
    DI.IterationStepFullProblem(opts);    
	
    opts = struct('noIterations',20,'lambda',0.8,'solveSquare',false);
    DI.IterationStepFullProblem(opts);    
    
	%opts = struct('noIterations',20,'lambda',0.8,'solveSquare',true);
    %DI.IterationStepFullProblem(opts);    
    
    DI.FindStagnationPoint();
    DI.SavePlotResults();
	DI.PlotErrorIterations();
                                     
    
end
