function SeppecherSingleFluid()

    close all;
    %% Parameters    
    PhysArea = struct('N',[40,30],'y2Min',0,'y2Max',20,'L1',10,'IntInterval',[-10,10]);

    PlotArea = struct('y1Min',-15,'y1Max',15,'N1',100,'N2',100);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('theta',pi/2,'g',0,...
                        'D_A',0,...
                        'zeta',10+2/3,...
                        'Cak',0.05,'Cn',4/3,...
                        'UWall',-1,...
                        'phi_m',3,...
                        'nParticles',0);
                    
    config = v2struct(optsPhys,optsNum);   
            
    DI = DiffuseInterfaceSingleFluid(config);
    DI.Preprocess();                  
    DI.IterationStepFullProblem_Seppecher();
                             
    
end
