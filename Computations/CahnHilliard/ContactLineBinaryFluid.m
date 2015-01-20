function ContactLineBinaryFluid
    close all;
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_InnerRegion'],'ORG');    

    PhysArea = struct('N',[50,40],'y2Min',0,'y2Max',20,...
                      'L1',7,'IntInterval',[-10,10]);%,'NBorder',[30,200,30,200]);

    PlotArea = struct('y1Min',-15,'y1Max',15,'N1',80,'N2',80);   
    SubArea  = struct('shape','Box','N',[60,60],...
                      'y1Min',-5,'y1Max',5,'y2Min',0,'y2Max',10);   	

    optsNum  = v2struct(PhysArea,PlotArea);                   	

    optsPhys = struct('thetaEq',pi/2,...       
                       'theta',90*pi/180,...
                       'Cak',0.01,'Cn',1,...
                       'UWall',1,...                       
                       'mobility',10,...
                       'nParticles',0);

    config = v2struct(optsPhys,optsNum);   

    opts = struct('noIterations',20,'lambda',0.8,'Seppecher_red',1);


    pars.Cak   = (0.01:0.01:0.03)';
    pars.y2Max = (10:5:15);

    dataM = DataStorage('NumericalExperiment',@RunNumericalExperiment,pars,config);

    function dataM = RunNumericalExperiment(pars,config)

        for i = 1:length(pars.Cak)
            config.optsPhys.Cak           = pars.Cak(i);        
            for j = 1:length(pars.y2Max)            
                config.optsNum.PhysArea.y2Max = pars.y2Max(j);

                DI = DiffuseInterfaceBinaryFluid(config);
                DI.Preprocess();

                DI.IterationStepFullProblem(opts);    

                opts.Seppecher_red = 2;
                opts.lambda        = 0.6;    
                %DI.optsPhys.mobility = m;
                DI.IterationStepFullProblem(opts);    
                DI.PlotResults();	  


                dataM(i,j).config = config;                
                dataM(i,j).theta                   = DI.GetThetaY2();
                dataM(i,j).y2                      = DI.IDC.Pts.y2;
                [dataM(i,j).lambda,dataM(i,j).Cin] = DI.FitSliplength();

                clear('DI');
            end
        end
    end

    thetaEq = config.optsPhys.thetaEq;

    for i = 1:size(dataM,1)

        Ca     = 3/4*pars.Cak(i);
        figure('name',['Ca = ',num2str(Ca)]);            
            
        for jn = 1:size(dataM,2)
            
            L      = pars.y2Max(j);
            y2     = dataM(i,j).y2;      
            y2P    = y2(2:end);
            
            lambda = (dataM(i,j).lambda);
            Cin    = dataM(i,j).Cin;            
            theta_Ana = GHR_Inv(Ca*log(y2P/lambda)+GHR_lambdaEta(thetaEq,1),1) + Ca*Cin;
            plot(y2,180/pi*dataM(i,j).theta,'ok'); hold on;
            plot(y2P,180/pi*theta_Ana,'m','linewidth',1.3); hold on;
                                     
            disp(['lambda = ',num2str(lambda),...
                  ' , Cin = ',num2str(Cin)]);
        end
    end

    % 
    % opts.solveSquare   = false;
    % DI.IterationStepFullProblem(opts);    

    % DI.optsNum.SubArea = SubArea;
    % DI.Preprocess_SubArea();
    % DI.PostProcess_Flux;
    % DI.IDC.SetUpBorders([30,1000,30,200]);
    % DI.FindAB();
    % 
    % %DI.FindStagnationPoint();
    % DI.PlotResults();	               
end