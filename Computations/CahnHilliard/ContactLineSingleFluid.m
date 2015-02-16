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
    pars.Cak    = (0.01:0.005:0.015)';
    pars.y2Max  = (16:2:24);             
    dataM       = DataStorage('NumericalExperiment',@RunNumericalExperiment,pars,[]);
    
    
    
    cols = {'g','b','c','k','r'};  nocols = length(cols);
	syms = {'<','d','s','o','>'};  nosyms = length(syms);
	
    
    PlotAsymptoticResults(dataM,'hatL','\hat{L}');
    PlotAsymptoticResults(dataM,'stagnationPointY2','y_{S2}');
    
    CompareAllData(dataM,1);
    CompareAllData(dataM,2);
            
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
                DI.ComputeInterfaceContour();
               % DI.PlotResults();	  
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
                                     
       
    function PlotAsymptoticResults(dataM,parameter,parName)
        figure('Position',[0 0 800 800],'color','white');
         
        for j1 = 1:size(dataM,1)
            for j2 = 1:size(dataM,2)
                y2Max(j2)    = dataM(j1,j2).config.optsNum.PhysArea.y2Max;
                par(j2)      = dataM(j1,j2).(parameter);                
            end                
            col = cols{mod(j2,nocols)+1};
            plot(y2Max,par,...
                ['-',syms{mod(j1,nosyms)+1},col],...
                'MarkerSize',8,'MarkerFaceColor',col); hold on;                
        end
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y_{2,max}$','Interpreter','Latex','fontsize',20);
        ylabel(['$',parName,'$'],'Interpreter','Latex','fontsize',20);        
        %ylim([90,100]);        
        SaveFigure([parameter,'_vs_l_diff'],pars);
    end
    function CompareAllData(dataM,j1)
        dataM = {dataM};
        thetaEq = dataM{1}(1,1).config.optsPhys.thetaEq;        

        figure('Position',[0 0 800 800],'color','white');
        hatL_M = zeros(size(dataM));
        Sy2_M  = zeros(size(dataM));
        Ca     = 3/4*pars.Cak(j1);
        
        for j0 = 1:length(dataM)
            for j2 = 1:size(dataM{j0},2)
                hatL_M(j1,j2) = dataM{j0}(j1,j2).hatL;         
                Sy2_M(j1,j2)  = dataM{j0}(j1,j2).stagnationPointY2;         
                PlotInterfaceSlope(dataM{j0},j1,j2,cols{mod(j0,nocols)+1},syms{mod(j2,nosyms)+1});
            end                                       
        end       
        

        y2P       = (1:0.1:max(dataM{end}(j1,j2).y2 - 4))';
        hatL      = dataM{end}(j1,j2).hatL;
        %theta_Ana = GHR_Inv(Ca*log(y2P/hatL)+GHR_lambdaEta(thetaEq,1),1);
        %plot(y2P,180/pi*theta_Ana,'m','linewidth',3); hold on;        


        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta[^\circ]$','Interpreter','Latex','fontsize',20);        
        %ylim([90,100]);
        
        SaveFigure(['InterfaceSlope_Ca_',num2str(Ca)],pars);
        
    end

    function PlotInterfaceSlope(dataM,i1,i2,col,sym)                    
        %hatL   = 0.46*l_diff;%46
        y2     = dataM(i1,i2).y2;
        mark   = (y2 < (max(y2)-4));
        plot(y2(mark),180/pi*dataM(i1,i2).theta(mark),['-',col,sym],'MarkerSize',5,'MarkerFaceColor',col); hold on;                              
    end
    
end
