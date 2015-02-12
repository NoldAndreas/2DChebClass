function ContactLineSingleFluid()

    close all;
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'DI_SingleFluid_InnerRegion'],'ORG');    
        
    PhysArea = struct('N',[50,40],'y2Min',0,'y2Max',15,...
                      'L1',7,'IntInterval',[-10,10]);%,'NBorder',[30,200,30,200]);

    PlotArea = struct('y1Min',-15,'y1Max',15,'N1',80,'N2',80);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('thetaEq',pi/2,...  
                       'zeta',10,...
                       'Cak',0.02,'Cn',1,...
                       'UWall',1,...
                       'phi_m',4,...
                       'nParticles',0);
                    
    
    pars.config = v2struct(optsPhys,optsNum);
    pars.Cak    = (0.01:0.005:0.02)';
    pars.y2Max  = (16:2:24);             
    dataM = DataStorage('NumericalExperiment',@RunNumericalExperiment,pars,[]);
        
            
%     DI = DiffuseInterfaceSingleFluid(config);
%     DI.Preprocess();                  
%     DI.IterationStepFullProblem(30);    
%     DI.FindStagnationPoint();    
% 	DI.PlotErrorIterations();
    
     function dataM = RunNumericalExperiment(pars,h)

        config = pars.config;
        for i = 1:length(pars.Cak)
            config.optsPhys.Cak           = pars.Cak(i);
            for j = 1:length(pars.y2Max)            
                config.optsNum.PhysArea.y2Max = pars.y2Max(j);

                DI = DiffuseInterfaceSingleFluid(config);
                DI.Preprocess();               
                DI.IterationStepFullProblem(40);                                
                DI.PlotResults();	  
                DI.FindStagnationPoint();

                dataM(i,j).config   = config;                
                dataM(i,j).Ca       = 3/4*pars.Cak(i);
                dataM(i,j).theta    = DI.GetThetaY2();
                dataM(i,j).y2       = DI.IDC.Pts.y2;
                dataM(i,j).hatL     = DI.FitSliplength();
                dataM(i,j).stagnationPointY2 = DI.StagnationPoint.y2_kv(1);

                close all;
                clear('DI');
            end
        end
    end
                                     
    
end
