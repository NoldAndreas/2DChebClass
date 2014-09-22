function ShearWithInterfaceBinaryFluid()

    close all;
    
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_ShearedInterface'],'ORG');
    %% Parameters    
    PhysArea = struct('N',[40,30],'y2Min',0,'y2Max',10,'L1',5,... %12,80,50
                      'NBorder',200);

    PlotArea = struct('y1Min',-20,'y1Max',20,'N1',80,...
                      'y2Min',0,'y2Max',PhysArea.y2Max,'N2',80);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('thetaEq',[pi/2,pi/2],...   
                      'eta',1,...
                      'Cak',0.05,'Cn',1,...
                      'UWall',[1,-1],...                      .
                      'mobility',10,...
                      'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
                      
    for y2M = 10%30:5:40%15:2.5:30
        
        config.optsNum.PhysArea.y2Max = y2M;
        config.optsNum.PlotArea.y2Max = y2M;
        
        DI = DiffuseInterfaceBinaryFluid(config);
        DI.Preprocess();                             
        DI.SolveMovingContactLine(60);    
        DI.SavePlotResults();
        DI.PlotErrorIterations();
        clear 'DI'
    end


end
