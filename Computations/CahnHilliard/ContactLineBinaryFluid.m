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

    pars.Cak   = (0.005:0.005:0.01)';
    pars.y2Max = (12:2:18);

    dataM = DataStorage('NumericalExperiment',@RunNumericalExperiment,pars,config);

    function dataM = RunNumericalExperiment(pars,config)

        for i = 1:length(pars.Cak)
            config.optsPhys.Cak           = pars.Cak(i);        
            for j = 1:length(pars.y2Max)            
                config.optsNum.PhysArea.y2Max = pars.y2Max(j);

                DI = DiffuseInterfaceBinaryFluid(config);
                DI.Preprocess();

                DI.IterationStepFullProblem(opts);    
                DI.GetYueParameters(); 

                opts.Seppecher_red = 2;
                opts.lambda        = 0.6;    
                %DI.optsPhys.mobility = m;
                DI.IterationStepFullProblem(opts);    
                DI.PlotResults();	  

                dataM(i,j).config                  = config;                
                dataM(i,j).theta                   = DI.GetThetaY2();
                dataM(i,j).y2                      = DI.IDC.Pts.y2;
                [dataM(i,j).lambda,dataM(i,j).Cin] = DI.FitSliplength();

                close all;
                clear('DI');
            end
        end
    end

    thetaEq = config.optsPhys.thetaEq;
    cols = {'g','b','m','k','r'};
    syms = {'d','s','o','>','<'};

    for i1 = 1:size(dataM,1)

        Ca     = 3/4*pars.Cak(i1);
        figure('name',['Ca = ',num2str(Ca)],...
                   'Position',[0 0 800 800],...
                   'color','white');            
            
        for i2 = 1:size(dataM,2)
            
            y2     = dataM(i1,i2).y2;      
            L      = max(y2);
            y2P    = y2(2:end);
            
            lambda = (dataM(i1,i2).lambda);
            Cin    = dataM(i1,i2).Cin;            
            theta_Ana = GHR_Inv(Ca*log(y2P/lambda)+GHR_lambdaEta(thetaEq,1),1) + Ca*Cin;
            plot(y2,180/pi*dataM(i1,i2).theta,['k',syms{i2}],'MarkerSize',7,'MarkerFaceColor','k'); hold on;
            plot(y2P,180/pi*theta_Ana,'m','linewidth',1.5); hold on;
                                     
            disp(['lambda = ',num2str(lambda),...
                  ' , Cin = ',num2str(Cin)]);
        end
        
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',15);
        xlabel('$y_2$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta$','Interpreter','Latex','fontsize',20);        

        filename = ['InterfaceSlope_Ca_',num2str(Ca)];
        print2eps([dirData filesep filename],gcf);
        saveas(gcf,[dirData filesep filename '.fig']);        
        disp(['Figures saved in ',dirData filesep filename '.fig/eps']);

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