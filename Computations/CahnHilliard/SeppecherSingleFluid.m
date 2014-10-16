function SeppecherSingleFluid()

    close all;
    %% Parameters    
    PhysArea = struct('N',[50,30],'y2Min',0,'y2Max',15,'L1',10,'IntInterval',[-10,10]);

    PlotArea = struct('y1Min',-15,'y1Max',15,'N1',80,'N2',80);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('thetaEq',pi/2,...  
                       'zeta',10,...
                       'Cak',0.03,'Cn',1,...
                       'UWall',1,...
                       'phi_m',4,...
                       'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
            
    DI = DiffuseInterfaceSingleFluid(config);
    DI.Preprocess();                  
    DI.IterationStepFullProblem_Seppecher(40);
    DI.SavePlotResults();
	DI.PlotErrorIterations();
                                     
    
end
