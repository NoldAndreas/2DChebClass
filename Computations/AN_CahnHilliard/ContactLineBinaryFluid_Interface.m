    close all;    
    AddPaths('DIBinaryPaper');            

    PhysArea = struct(...%'shape','InfCapillary',...%_Interface',...
                      'N',[50,40],'y2Min',0,'y2Max',40,...
                      'L1',7,'IntInterval',[-10,10]);%,'NBorder',[30,200,30,200]);

    PlotArea = struct('y1Min',-5,'y1Max',5,'N1',80,'N2',80);   
    SubArea  = struct('shape','Box','N',[60,60],...
                      'y1Min',-5,'y1Max',5,'y2Min',0,'y2Max',10);   	
                  
    %optsSolv = struct('noIterations',40,'lambda',0.8,'Seppecher_red',1);

    optsNum  = v2struct(PhysArea,PlotArea);                   	

    optsPhys = struct('thetaEq',pi/2,...       
                       'theta',90*pi/180,...
                       'Cak',0.01,'Cn',1,...
                       'UWall',1,...                       
                       'l_diff',2.0,...%'mobility',10,...
                       'nParticles',0);
                                          
    config = v2struct(optsPhys,optsNum);

    DI = DiffuseInterfaceBinaryFluid(config);
    DI.Preprocess();

	DI.IterationStepFullProblem();                    
	DI.PostProcess(); 