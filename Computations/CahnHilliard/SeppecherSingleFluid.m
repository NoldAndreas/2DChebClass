function SeppecherSingleFluid()

    close all;
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_InnerRegion'],'ORG');    
        
    PhysArea = struct('N',[50,30],'y2Min',0,'y2Max',10,...
                      'L1',7,'IntInterval',[-10,10]);%,'NBorder',[30,200,30,200]);

    PlotArea = struct('y1Min',-15,'y1Max',15,'N1',80,'N2',80);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('thetaEq',pi/2,...  
                       'zeta',10,...
                       'Cak',0.02,'Cn',1,...
                       'UWall',1,...
                       'phi_m',4,...
                       'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
            
    DI = DiffuseInterfaceSingleFluid(config);
    DI.Preprocess();                  
    DI.IterationStepFullProblem_Seppecher(20);
    DI.SavePlotResults();
	DI.PlotErrorIterations();
                                     
    
end
