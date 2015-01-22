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
                       'l_diff',0.5,...%'mobility',10,...
                       'nParticles',0);

    config = v2struct(optsPhys,optsNum);   

    opts = struct('noIterations',20,'lambda',0.8,'Seppecher_red',1);

    pars.config = config;
    pars.Cak   = (0.005:0.0025:0.015)';
    pars.y2Max = (12:2:18);

    dataM = DataStorage('NumericalExperiment',@RunNumericalExperiment,pars,[]);

    function dataM = RunNumericalExperiment(pars,h)

        config = pars.config;
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
                
                DI.GetYueParameters(); 
                DI.PlotResults();	  

                dataM(i,j).config   = config;                
                dataM(i,j).theta    = DI.GetThetaY2();
                dataM(i,j).y2       = DI.IDC.Pts.y2;
                dataM(i,j).hatL     = DI.FitSliplength();
                dataM(i,j).stagnationPointY2 = DI.StagnationPoint.y2_kv(1);

                close all;
                clear('DI');
            end
        end
    end

    thetaEq = config.optsPhys.thetaEq;
    cols = {'b','m','k','r','g'};
    syms = {'d','s','o','>','<'};
    
    figure('Position',[0 0 800 800],'color','white');
    
    for i1 = 1:size(dataM,1)
        
        Ca          = 3/4*pars.Cak(i1);   
        hatL_Av(i1) = 0;
        
        for i2 = 1:size(dataM,2)
            hatL_Av(i1) = hatL_Av(i1) + (dataM(i1,i2).hatL);
            PlotInterfaceSlope(i1,i2,cols{i1},syms{i2});
        end        
        hatL_Av(i1) = hatL_Av(i1)/size(dataM,2);
        
        y2P       = dataM(i1,i2).y2(2:end);
        theta_Ana = GHR_Inv(Ca*log(y2P/hatL_Av(i1))+GHR_lambdaEta(thetaEq,1),1);
        plot(y2P,180/pi*theta_Ana,cols{i1},'linewidth',1.5); hold on;        

    end
    
	set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);
    xlabel('$y$','Interpreter','Latex','fontsize',20);
    ylabel('$\theta[^\circ]$','Interpreter','Latex','fontsize',20);        
    ylim([90,105]);
    
	filename = ['InterfaceSlope'];
    print2eps([dirData filesep filename],gcf);
    saveas(gcf,[dirData filesep filename '.fig']);        
    disp(['Figures saved in ',dirData filesep filename '.fig/eps']);

    figure('Position',[0 0 800 800],'color','white');
    plot(pars.Cak,hatL_Av.^2,'o-','MarkerSize',7,'MarkerFaceColor','k','linewidth',1.5);
    xlabel('Ca','Interpreter','Latex','fontsize',20);       
    ylabel('$(\hat L)^2$','Interpreter','Latex','fontsize',20);       

	set(gca,'linewidth',1.5);
	set(gca,'fontsize',15);
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
    
    function PlotInterfaceSlope(i1,i2,col,sym)
                    
        y2     = dataM(i1,i2).y2;                          
        hatL   = (dataM(i1,i2).hatL);                        
        
        plot(y2,180/pi*dataM(i1,i2).theta,[col,sym],'MarkerSize',7,'MarkerFaceColor',col); hold on;                
        disp(['hatL = ',num2str(hatL)]);                
    end
end