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
                  
    %optsSolv = struct('noIterations',40,'lambda',0.8,'Seppecher_red',1);

    optsNum  = v2struct(PhysArea,PlotArea);                   	

    optsPhys = struct('thetaEq',pi/2,...       
                       'theta',90*pi/180,...
                       'Cak',0.01,'Cn',1,...
                       'UWall',1,...                       
                       'l_diff',0.5,...%'mobility',10,...
                       'nParticles',0);

    parameters.config = v2struct(optsPhys,optsNum);
    parameters.Cak   = [0.005;0.01];%(0.005:0.0025:0.01)';
    parameters.y2Max = 24;%(16:2:24);            
    parameters.l_d   = 1:0.25:3.0;%0.25:0.25:2.5;        
%         for k = 1:length(l_d)
%             parameters.config.optsPhys.l_diff = params.l_d(k);    
%             parameters.y2Max                  = y2Max*parameters.config.optsPhys.l_diff;
%             
%         end
	[dataM,~,res] = DataStorage('NumericalExperiment',@RunNumericalExperiment,parameters,[]);
    dataN = Rescale(dataM);clear('dataM');
    
    
    cols = {'g','b','c','k','r'};  nocols = length(cols);
	syms = {'<','d','s','o','>'};  nosyms = length(syms);            
    
    SaveFigures(dataN,res,parameters);
    
 %   PlotAsymptoticResults_Y2Max(dataM,'IsoInterface.hatL','\hat{L}');
%     PlotAsymptoticResults_Y2Max(dataM,'stagnationPointY2','y_{S2}');
%     PlotAsymptoticResults_Y2Max(dataM,'muPlusInf','muPlusInf');
%     PlotAsymptoticResults_Y2Max(dataM,'muMinusInf','muMinusInf');
%     PlotAsymptoticResults_Y2Max(dataM,'muMaxAbs','muMaxAbs');
%     PlotAsymptoticResults_Y2Max(dataM,'pMin','pMin');
%     PlotAsymptoticResults_Y2Max(dataM,'pMax','pMax');

  %  PlotAsymptoticInterfaceResults(dataN,1,'theta',{'IsoInterface'},'\theta');    
    PlotAsymptoticInterfaceResults(dataN,1,'kappa',{'IsoInterface'},'\kappa');    
    PlotAsymptoticInterfaceResults(dataN,1,'mu',{'IsoInterface'},'\mu');
    PlotAsymptoticInterfaceResults(dataN,1,'mu_ddy2',{'IsoInterface'},'\frac{\partial^2\mu}{\partial y_2^2}');
    PlotAsymptoticInterfaceResults(dataN,1,'p',{'IsoInterface'},'p');
    PlotAsymptoticInterfaceResults(dataN,1,'kappa_Cakmu',{'IsoInterface'},'\kappa/(Ca_k\mu)');    
    
    
	PlotAsymptoticInterfaceResults(dataN,1,'mu',{'mu_y2'},'\mu');        
    PlotAsymptoticInterfaceResults(dataN,1,'mu_ddy2',{'mu_y2'},'\frac{\partial^2\mu}{\partial y_2^2}');                
    
    PlotAsymptoticResults(dataN,'hatL',{'IsoInterface'}); 
    PlotAsymptoticResults(dataN,'stagnationPointY2');
    PlotAsymptoticResults(dataN,'muPlusInf');
    PlotAsymptoticResults(dataN,'muMinusInf');
    PlotAsymptoticResults(dataN,'muMaxAbs');
    PlotAsymptoticResults(dataN,'pMin');
    PlotAsymptoticResults(dataN,'pMax');
                    
   % CompareAllData(dataM,1);
   %CompareAllData(dataM,2);
   %CompareAllData(dataM,3);
    
    %dataM{l_diff}(Cak,y2Max)
    
    function SaveFigures(dataN,res,params) 
        params.comment = 'run by ContactLineBinaryFluid.m';
        [~,fn]   = fileparts(res.Filename);        
        filename = ['NumericalExperiment' filesep fn '_'];
        
        
        PlotAsymptoticInterfaceResults(dataN,1,'u_n',{'IsoInterface','noLegend'},'u_n');    
        plot([0,24],[0 0],'k--','linewidth',1.5);
        SaveFigure([filename, 'Interface_u_n'],params);        
        
        PlotAsymptoticInterfaceResults(dataN,1,'u_t',{'IsoInterface','noLegend'},'u_t');    
        plot([0,24],[0 0],'k--','linewidth',1.5);
        SaveFigure([filename, 'Interface_u_t'],params);        
        
        PlotAsymptoticInterfaceResults(dataN,1,'flux_n',{'IsoInterface','noLegend'},'j_n');    
        plot([0,24],[0 0],'k--','linewidth',1.5);
        SaveFigure([filename, 'Interface_flux_n'],params);        
        
        PlotAsymptoticInterfaceResults(dataN,1,'flux_t',{'IsoInterface','noLegend'},'j_t');    hold on;
        plot([0,24],[0 0],'k--','linewidth',1.5);
        SaveFigure([filename, 'Interface_flux_t'],params);                
        
        PlotAsymptoticInterfaceResults(dataN,1,'mu',{'mu_y2','noLegend'},'\mu');        
        SaveFigure([filename, 'wall_mu'],params);
        
        PlotAsymptoticInterfaceResults(dataN,1,'mu_dy1',{'mu_y2','noLegend'},'\frac{\partial\mu}{\partial y_1}');                
        SaveFigure([filename, 'wall_mu_dy1'],params);
        
        PlotAsymptoticInterfaceResults(dataN,1,'mu_ddy2',{'mu_y2','noLegend'},'\frac{\partial^2\mu}{\partial y_2^2}');                
        SaveFigure([filename, 'wall_mu_ddy2'],params);
        
        PlotAsymptoticResults(dataN,'hatL',{'IsoInterface'},'$\hat L$'); 
        SaveFigure([filename, 'hatL'],params);
                       
        PlotAsymptoticResults(dataN,'stagnationPointY2',{},'$y_{2,S}$');         
        SaveFigure([filename, 'stagnationPointY2'],params);
        
        PlotAsymptoticInterfaceResults(dataN,1,'kappa',{'IsoInterface','noLegend'},'\kappa');    
        SaveFigure([filename, 'Interfacekappa'],params);
                
        PlotAsymptoticInterfaceResults(dataN,1,'mu',{'IsoInterface','noLegend'},'\mu');
        SaveFigure([filename, 'InterfaceMu'],params);
        
        %PlotAsymptoticInterfaceResults(dataN,1,'mu_ddy2',{'IsoInterface'},'\frac{\partial^2\mu}{\partial y_2^2}');
        PlotAsymptoticInterfaceResults(dataN,1,'p',{'IsoInterface','noLegend'},'p');
        SaveFigure([filename, 'InterfaceP'],params);        
        
        PlotAsymptoticInterfaceResults(dataN,1,'kappa_Cakmu',{'IsoInterface','noLegend'},'\kappa/(Ca_k\mu)');     hold on;
        plot([0,24],[1.5 1.5],'k--','linewidth',1.5);
        SaveFigure([filename, 'Interface_KappaCak_Mu'],params);
    end
    function PlotAsymptoticInterfaceResults(dataM,i_Cak,value,opts,ylabelStr)
        if(nargin < 4)
            opts = {};                        
        end        
        if(ischar(opts))
            opts = {opts};
        end
        legendstr = {};
        figure('Position',[0 0 800 800],'color','white');
        title(['Ca = ',num2str(dataM{1}(i_Cak,1).Ca)]);
        for j2 = 1:size(dataM{1},2)
            for j0 = 1:length(dataM)                
                Cn    = dataM{j0}(i_Cak,j2).Cn;                
                col       = cols{mod(j0,nocols)+1};
                                
                if(IsOption(opts,'IsoInterface'))
                    val = dataM{j0}(i_Cak,j2).IsoInterface.(value); 
                    y  = dataM{j0}(i_Cak,j2).IsoInterface.y2;
                elseif(IsOption(opts,'mu_y2'))
                    val = dataM{j0}(i_Cak,j2).mu_y2{1}.(value);
                    y  = dataM{j0}(i_Cak,j2).mu_y2{1}.pts.y1_kv;                                        
                end
                
                if(IsOption(opts,'SC'))
                    plot(y,val,['-',syms{mod(j2,nosyms)+1},col],...
                        'MarkerSize',8,'MarkerFaceColor',col);
                else
                    plot(y,val,['-',col],'linewidth',1.5); 
                end
                hold on;                                
                if(size(dataM{1},2)>1)
                    legendstr(end+1) = {['Cn = ',num2str(Cn),' y_{2,Max} = ',num2str(dataM{j0}(i_Cak,j2).y2Max)]};
                else
                    legendstr(end+1) = {['Cn = ',num2str(Cn)]};
                end
            end
        end  
        if(~IsOption(opts,'noLegend'))            
            legend(legendstr);
        end
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);        
                       
        if(IsOption(opts,'mu_y2'))
            xlabel('$y_1$','Interpreter','Latex','fontsize',20);
            ylabel(['$',ylabelStr,'$ at $y_2 = ',num2str(dataM{1}(i_Cak,1).mu_y2{1}.y2),'$'],'Interpreter','Latex','fontsize',20);  
        else            
            xlabel('$y_2$','Interpreter','Latex','fontsize',20);
            ylabel(['$',ylabelStr,'$'],'Interpreter','Latex','fontsize',20);  
        end        
                
    end        
    function PlotAsymptoticResults(dataM,parameter,opts,ylabelStr)        
        
        if(nargin == 2)
            opts = {};                        
        end        
        if(ischar(opts))
            opts = {opts};
        end
        
        
        figure('Position',[0 0 800 800],'color','white');
        legendstr = {}; 
        for i_Cak = 1:size(dataM{1},1)
            Ca = dataM{1}(i_Cak,1).Ca;
            for j2 = 1:size(dataM{1},2)
                for j0 = 1:length(dataM)
                    Cn(j0)    = dataM{j0}(i_Cak,j2).Cn;
                    if(IsOption(opts,'IsoInterface'))
                        par(j0)      = dataM{j0}(i_Cak,j2).IsoInterface.(parameter);                
                    else                        
                        par(j0)      = dataM{j0}(i_Cak,j2).(parameter);                
                    end
                end 
                col = cols{mod(i_Cak,nocols)+1};                
                sym = syms{mod(j2,nosyms)+1};
                                
                plot(Cn,par,...
                    ['-',sym,col],...
                    'MarkerSize',8,'MarkerFaceColor',col); hold on;                
                
                legendstr(end+1) = {['Ca = ',num2str(Ca),' y_{2,Max} = ',num2str(dataM{1}(i_Cak,j2).y2Max)]};
             end  
        end
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        if(~IsOption(opts,'noLegend'))            
            legend(legendstr);
        end
        xlabel('Cn','Interpreter','Latex','fontsize',20);
        if(nargin >= 4)
            ylabel(ylabelStr,'Interpreter','Latex','fontsize',20);                
        else
            ylabel(parameter,'Interpreter','Latex','fontsize',20);                
        end
        %ylim([90,100]);        
        %SaveFigure([parameter,'_vs_l_diff']);
    end        
    function PlotAsymptoticResults_Y2Max(dataM,parameter,parName)
        figure('Position',[0 0 800 800],'color','white');
         
        for j1 = 1:size(dataM{1},1)
            for j0 = 1:length(dataM)
                l_diff = dataM{j0}(j1,1).config.optsPhys.l_diff;
                for j2 = 1:size(dataM{1},2)
                    y2M(j2)  = dataM{j0}(j1,j2).config.optsNum.PhysArea.y2Max;
                    par(j2)    = dataM{j0}(j1,j2).(parameter);                
                end                 
                col = cols{mod(j0,nocols)+1};
                plot(y2M/l_diff,par,...
                    ['-',syms{mod(j1,nosyms)+1},col],...
                    'MarkerSize',8,'MarkerFaceColor',col); hold on;                
             end  
        end
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y_{2,max}/l_d$','Interpreter','Latex','fontsize',20);
        ylabel(['$',parName,'/ \ell_D$'],'Interpreter','Latex','fontsize',20);        
        %ylim([90,100]);
        SaveFigure([parameter,'_vs_l_diff']);
    end
    function dataM = RunNumericalExperiment(pars,h)
    
        config = pars.config;
        for k0 = 1:length(pars.l_d)
            
            l_diff                 = pars.l_d(k0);
            config.optsPhys.l_diff = l_diff;    
            
            for j = 1:length(pars.y2Max)            

                config.optsNum.PhysArea.y2Max = pars.y2Max(j)*l_diff;

                DI = DiffuseInterfaceBinaryFluid(config);
                DI.Preprocess();
                for i = 1:length(pars.Cak)

                    DI.SetCak(pars.Cak(i));                
                    DI.IterationStepFullProblem();                    
                    DI.PostProcess();                       

                    dataM{k0}(i,j).config.optsNum    = DI.optsNum;
                    dataM{k0}(i,j).config.optsPhys   = DI.optsPhys;
                    dataM{k0}(i,j).Ca                = 3/4*pars.Cak(i);                                
                    dataM{k0}(i,j).muMinusInf        = DI.mu(DI.IDC.Ind.left & DI.IDC.Ind.bottom);
                    dataM{k0}(i,j).muPlusInf         = DI.mu(DI.IDC.Ind.right & DI.IDC.Ind.bottom);
                    dataM{k0}(i,j).muMaxAbs          = max(abs(DI.mu));
                    dataM{k0}(i,j).pMax              = max((DI.p));
                    dataM{k0}(i,j).pMin              = min((DI.p));                              
                    dataM{k0}(i,j).IsoInterface      = DI.IsoInterface;                                
                    dataM{k0}(i,j).stagnationPointY2 = DI.StagnationPoint.y2_kv(1);                    

                    y2 = 0:0.5:3;
                    for kk = 1:length(y2)
                        [mu,pts]      = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.mu);                    
                        [mu_dy1,pts]  = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy1*DI.mu);
                        [mu_ddy2,pts] = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.DDy2*DI.mu);
                        dataM{k0}(i,j).mu_y2{kk} = struct('mu',mu,...
                                                      'mu_dy1',mu_dy1,...
                                                      'mu_ddy2',mu_ddy2,...
                                                      'pts',pts,'y2',y2(kk));                        
                    end

                    close all;
                end
                clear('DI');
            end
        end        
    end  
    function CompareAllData(dataM,j1)
        thetaEq = dataM{1}(1,1).config.optsPhys.thetaEq;        

        figure('Position',[0 0 800 800],'color','white');
        hatL_M = zeros(size(dataM));
        Sy2_M  = zeros(size(dataM));
        Ca     = 3/4*dataM{1}(j1,1).config.optsPhys.Cak;
        
        for j0 = 1:length(dataM)
            for j2 = 1:size(dataM{j0},2)
                hatL_M(j1,j2) = dataM{j0}(j1,j2).hatL;         
                Sy2_M(j1,j2)  = dataM{j0}(j1,j2).stagnationPointY2;         
                PlotRescaledInterfaceSlope(dataM{j0},j1,j2,cols{mod(j0,nocols)+1},syms{mod(j2,nosyms)+1});
            end                                       
        end       
        
        l_diff    = dataM{end}(j1,j2).config.optsPhys.l_diff;
        y2P       = (1:0.1:max(dataM{end}(j1,j2).y2 - 4*l_diff))'/l_diff;
        hatL      = dataM{end}(j1,j2).hatL/l_diff;
        theta_Ana = GHR_Inv(Ca*log(y2P/hatL)+GHR_lambdaEta(thetaEq,1),1);
        plot(y2P,180/pi*theta_Ana,'m','linewidth',3); hold on;        


        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y/\ell_{d}$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta[^\circ]$','Interpreter','Latex','fontsize',20);        
        %ylim([90,100]);
        
        SaveFigure(['InterfaceSlope_Ca_',num2str(Ca)],pars);
        
    end
    function PlotData(dataM)
        thetaEq = dataM(1,1).config.optsPhys.thetaEq;
        nocols = 5;
        nosyms = 5;
        cols = {'g','b','m','k','r'};
        syms = {'<','d','s','o','>'};

        figure('Position',[0 0 800 800],'color','white');
        hatL_M = zeros(size(dataM));
        Sy2_M  = zeros(size(dataM));

        for i1 = 1:size(dataM,1)

            Ca          = 3/4*pars.Cak(i1);           

            for i2 = 1:size(dataM,2)
                hatL_M(i1,i2) = dataM(i1,i2).hatL;         
                Sy2_M(i1,i2)  = dataM(i1,i2).stagnationPointY2;         
                PlotInterfaceSlope(dataM,i1,i2,cols{mod(i1,nocols)+1},syms{mod(i2,nosyms)+1});
            end        
            hatL_Av = sum(hatL_M(i1,:)/size(dataM,2));

            y2P       = dataM(i1,i2).y2(2:end);
            theta_Ana = GHR_Inv(Ca*log(y2P/hatL_Av)+GHR_lambdaEta(thetaEq,1),1);
            plot(y2P,180/pi*theta_Ana,cols{mod(i1,nocols)+1},'linewidth',1.5); hold on;        

        end       

        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta[^\circ]$','Interpreter','Latex','fontsize',20);        
        ylim([90,100]);
        
        
        l_diff = pars.config.optsPhys.l_diff;

        SaveFigure(['InterfaceSlope_l_d_',num2str(l_diff)],pars);
        %print2eps([dirData filesep filename],gcf);
        %saveas(gcf,[dirData filesep filename '.fig']);        
        %disp(['Figures saved in ',dirData filesep filename '.fig/eps']);

        hatL_Av = mean2(hatL_M); 
        Lb =  hatL_Av - min(min(hatL_M));
        Ub =  max(max(hatL_M))-hatL_Av;       
        disp(['hatL/l_diff = ',num2str(hatL_Av/l_diff),' +/- ',num2str(max(Lb,Ub)/l_diff)]);


        av = mean2(Sy2_M); 
        Lb =  av - min(min(Sy2_M));
        Ub =  max(max(Sy2_M))-av;           
        disp(['stagnation point y2/l_diff = ',num2str(av/l_diff),' +/- ',num2str(max(Lb,Ub)/l_diff)]);
    end
    
%     plot(pars.Cak,hatL_Av,'o-','MarkerSize',7,'MarkerFaceColor','k','linewidth',1.5);
%     errorbar(pars.Cak,hatL_Av,min(hatL_M,[],2)-hatL_Av,max(hatL_M,[],2)-hatL_Av);
%     %plot(pars.Cak,min(hatL_M,[],2),'o-','MarkerSize',7,'MarkerFaceColor','k','linewidth',1.5);
%     %plot(pars.Cak,max(hatL_M,[],2),'o-','MarkerSize',7,'MarkerFaceColor','k','linewidth',1.5);
%     xlabel('Ca','Interpreter','Latex','fontsize',20);       
%     ylabel('$\hat L$','Interpreter','Latex','fontsize',20);       
% 
% 	set(gca,'linewidth',1.5);
% 	set(gca,'fontsize',15);
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
    
    function PlotRescaledInterfaceSlope(dataM,i1,i2,col,sym)                    
        l_diff = dataM(i1,i2).config.optsPhys.l_diff;
        %hatL   = 0.46*l_diff;%46
        y2     = dataM(i1,i2).y2 /l_diff;
        mark   = (y2 < (max(y2)-4));
        plot(y2(mark),180/pi*dataM(i1,i2).theta(mark),[col,sym],'MarkerSize',5,'MarkerFaceColor',col); hold on;                              
    end   
    function PlotInterfaceSlope(dataM,i1,i2,col,sym)
                    
        y2     = dataM(i1,i2).y2;                          
        hatL   = (dataM(i1,i2).hatL);                        
        
        plot(y2,180/pi*dataM(i1,i2).theta,[col,sym],'MarkerSize',7,'MarkerFaceColor',col); hold on;                
        disp(['hatL = ',num2str(hatL)]);                
    end

    function dataN = Rescale(dataM)
        for k0 = 1:length(dataM)
            for k1 = 1:size(dataM{1},1)
                for k2 = 1:size(dataM{1},2)
                    Cn    = 1/dataM{k0}(k1,k2).config.optsPhys.l_diff;
                    Cak   = dataM{k0}(k1,k2).config.optsPhys.Cak;
                    
                    dataN{k0}(k1,k2).Cn  = Cn;           
                    dataN{k0}(k1,k2).Cak = Cak;
                    
                    dataN{k0}(k1,k2).y2Max  = dataM{k0}(k1,k2).config.optsNum.PhysArea.y2Max*Cn;

                    dataN{k0}(k1,k2).Ca              = dataM{k0}(k1,k2).Ca;                                
                    
                    
                    dataN{k0}(k1,k2).muMinusInf  = dataM{k0}(k1,k2).muMinusInf/(Cn*Cak);
                    dataN{k0}(k1,k2).muPlusInf   = dataM{k0}(k1,k2).muPlusInf/(Cn*Cak);
                    dataN{k0}(k1,k2).muMaxAbs    = dataM{k0}(k1,k2).muMaxAbs/(Cn*Cak);
                    dataN{k0}(k1,k2).pMax        = dataM{k0}(k1,k2).pMax/Cn;
                    dataN{k0}(k1,k2).pMin        = dataM{k0}(k1,k2).pMin/Cn;                                      
                    dataN{k0}(k1,k2).stagnationPointY2 = dataM{k0}(k1,k2).stagnationPointY2*Cn;
                    
                    for kk = 1:length(dataM{1}(1,1).mu_y2)
                        dataN{k0}(k1,k2).mu_y2{kk}.mu        = dataM{k0}(k1,k2).mu_y2{kk}.mu/(Cn*Cak);
                        dataN{k0}(k1,k2).mu_y2{kk}.mu_dy1    = dataM{k0}(k1,k2).mu_y2{kk}.mu_dy1/(Cn^2*Cak);
                        dataN{k0}(k1,k2).mu_y2{kk}.mu_ddy2   = dataM{k0}(k1,k2).mu_y2{kk}.mu_ddy2/(Cn^3*Cak);
                        dataN{k0}(k1,k2).mu_y2{kk}.pts.y1_kv = dataM{k0}(k1,k2).mu_y2{kk}.pts.y1_kv*Cn;
                        dataN{k0}(k1,k2).mu_y2{kk}.pts.y2_kv = dataM{k0}(k1,k2).mu_y2{kk}.pts.y2_kv*Cn;
                        dataN{k0}(k1,k2).mu_y2{kk}.y2        = dataM{k0}(k1,k2).mu_y2{kk}.y2*Cn;                        
                    end

                    
                    dataN{k0}(k1,k2).IsoInterface.y2       = dataM{k0}(k1,k2).IsoInterface.y2*Cn;
                    dataN{k0}(k1,k2).IsoInterface.theta    = dataM{k0}(k1,k2).IsoInterface.theta;
                    dataN{k0}(k1,k2).IsoInterface.hatL     = dataM{k0}(k1,k2).IsoInterface.hatL*Cn;
                    dataN{k0}(k1,k2).IsoInterface.mu_ddy2  = dataM{k0}(k1,k2).IsoInterface.mu_ddy2/(Cn^3*Cak);
                    dataN{k0}(k1,k2).IsoInterface.u_n  = dataM{k0}(k1,k2).IsoInterface.u_n;
                    dataN{k0}(k1,k2).IsoInterface.u_t  = dataM{k0}(k1,k2).IsoInterface.u_t;
                    dataN{k0}(k1,k2).IsoInterface.flux_n  = dataM{k0}(k1,k2).IsoInterface.flux_n;
                    dataN{k0}(k1,k2).IsoInterface.flux_t  = dataM{k0}(k1,k2).IsoInterface.flux_t;
                    
                    
                    dataN{k0}(k1,k2).IsoInterface.kappa = dataM{k0}(k1,k2).IsoInterface.kappa/Cn;
                    dataN{k0}(k1,k2).IsoInterface.p     = dataM{k0}(k1,k2).IsoInterface.p/Cn;
                    dataN{k0}(k1,k2).IsoInterface.mu    = dataM{k0}(k1,k2).IsoInterface.mu/(Cak*Cn);
                    dataN{k0}(k1,k2).IsoInterface.kappa_Cakmu       = dataM{k0}(k1,k2).IsoInterface.kappa./dataM{k0}(k1,k2).IsoInterface.mu;
                end
            end
        end
	end
end

