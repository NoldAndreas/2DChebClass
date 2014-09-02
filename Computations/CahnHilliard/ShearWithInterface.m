function ShearWithInterface()

    close all;
    
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_InnerRegion'],'ORG');    
    %% Parameters    
    PhysArea = struct('N',[70,40],'y2Min',0,'y2Max',25,'L1',10,... %12,80,50
                      'NBorder',200);

    PlotArea = struct('y1Min',-20,'y1Max',20,'N1',80,...
                      'y2Min',0,'y2Max',PhysArea.y2Max,'N2',80);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('thetaEq',[pi/2,pi/2],...   
                      'zeta',10,'eta',1,...
                      'Cak',0.1,'Cn',1,...
                      'UWall',[1,-1],...
                      'rho_m',4,...
                      'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
                      
    for y2M = 10:2.5:22.5
        config.PhysArea.y2Max = y2M;
        DI = DiffuseInterface(config);
        DI.Preprocess();                             
        DI.SolveMovingContactLine(15);    
        DI.SavePlotResults();
        DI.PlotErrorIterations();
        clear 'DI'
    end


end
