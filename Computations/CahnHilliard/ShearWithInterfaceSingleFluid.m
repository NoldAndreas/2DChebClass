function ShearWithInterfaceSingleFluid()

    close all;
    
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_ShearedInterface'],'ORG');
    %% Parameters    
    PhysArea = struct('N',[40,30],'y2Min',0,'y2Max',[],'L1',10,...
                      'IntInterval',[-10,10]);%,... %12,80,50
    PlotArea = struct('y1Min',-10,'y1Max',10,'N1',80,'N2',80);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('thetaEq',[pi/2,pi/2],...   
                      'zeta',10,...
                      'Cak',0.03,'Cn',1,...
                      'UWall',[1,-1],...
                      'phi_m',3,...
                      'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
                      
    for y2M = 15%30:5:40%15:2.5:30
        
        config.optsNum.PhysArea.y2Max = y2M;
        config.optsNum.PlotArea.y2Max = y2M;
        
        DI = DiffuseInterfaceSingleFluid(config);
        DI.Preprocess();         
        DI.IterationStepFullProblem(20);
        %DI.SolveMovingContactLine(20);    
        DI.SavePlotResults();
        DI.PlotErrorIterations();
        clear 'DI'
    end


end
