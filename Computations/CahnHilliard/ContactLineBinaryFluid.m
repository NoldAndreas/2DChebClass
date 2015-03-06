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
    parameters.Cak   = 0.005;%[0.005;0.01];%(0.005:0.0025:0.01)';
    y2Max            = 24;%(16:2:24);            
    l_d              = 3.0;%1:0.25:3.0;%0.25:0.25:2.5;
    
    for k = 1:length(l_d)
        parameters.config.optsPhys.l_diff = l_d(k);    
        parameters.y2Max                  = y2Max*parameters.config.optsPhys.l_diff;
        dataM{k} = DataStorage('NumericalExperiment',@RunNumericalExperiment,parameters,[],true);
    end
    
    for k0 = 1:length(dataM)
        for k1 = 1:size(dataM{1},1)
            for k2 = 1:size(dataM{1},2)
                
                l_d    = dataM{k0}(k1,k2).config.optsPhys.l_diff;     
                Cak    = dataM{k0}(k1,k2).config.optsPhys.Cak;     
                
                for k3 = 1:length(dataM{k0}(k1,k2).mu_y2)
                    Diff   = barychebdiff(dataM{k0}(k1,k2).mu_y2{k3}.pts.y1_kv/l_d,1);
                    dataM{k0}(k1,k2).mu_y2{k3}.mu_0P = Diff.Dx*dataM{k0}(k1,k2).mu_y2{k3}.mu;
                end
                
                dataM{k0}(k1,k2).IsoInterface.kappa_rescaled = dataM{k0}(k1,k2).IsoInterface.kappa*l_d;
                dataM{k0}(k1,k2).IsoInterface.p_rescaled     = dataM{k0}(k1,k2).IsoInterface.p*l_d;
                dataM{k0}(k1,k2).IsoInterface.mu_rescaled    = dataM{k0}(k1,k2).IsoInterface.mu*l_d;
                dataM{k0}(k1,k2).IsoInterface.kappa_mu       = dataM{k0}(k1,k2).IsoInterface.kappa./dataM{k0}(k1,k2).IsoInterface.mu;
                %dataM{k0}(k1,k2).IsoInterface.kappaRescaled = dataM{k0}(k1,k2).IsoInterface.kappa*l_d;
            end
        end
    end
                
    
	cols = {'g','b','c','k','r'};  nocols = length(cols);
	syms = {'<','d','s','o','>'};  nosyms = length(syms);
            
 %   PlotAsymptoticResults_Y2Max(dataM,'IsoInterface.hatL','\hat{L}');
%     PlotAsymptoticResults_Y2Max(dataM,'stagnationPointY2','y_{S2}');
%     PlotAsymptoticResults_Y2Max(dataM,'muPlusInf','muPlusInf');
%     PlotAsymptoticResults_Y2Max(dataM,'muMinusInf','muMinusInf');
%     PlotAsymptoticResults_Y2Max(dataM,'muMaxAbs','muMaxAbs');
%     PlotAsymptoticResults_Y2Max(dataM,'pMin','pMin');
%     PlotAsymptoticResults_Y2Max(dataM,'pMax','pMax');
    
    for ii = 1:length(dataM{1}(1,1).mu_y2)
        PlotAsymptoticInterfaceResults(dataM,1,ii,{'mu_y2'});    
    end
    
    PlotAsymptoticInterfaceResults(dataM,1,'kappa_rescaled',{'IsoInterface'});    
    PlotAsymptoticInterfaceResults(dataM,1,'mu_rescaled',{'IsoInterface'});
    PlotAsymptoticInterfaceResults(dataM,1,'p_rescaled',{'IsoInterface'});
    PlotAsymptoticInterfaceResults(dataM,1,'kappa_mu',{'IsoInterface'});
    
    
    PlotAsymptoticResults(dataM,'hatL',{'rescale','IsoInterface'});
    PlotAsymptoticResults(dataM,'stagnationPointY2','rescale');
    PlotAsymptoticResults(dataM,'muPlusInf',{'rescaleInv'});
    PlotAsymptoticResults(dataM,'muMinusInf',{'rescaleInv'});
    PlotAsymptoticResults(dataM,'muMaxAbs',{'rescaleInv','rescaleCa'});
    PlotAsymptoticResults(dataM,'pMin',{'rescaleInv'});
    PlotAsymptoticResults(dataM,'pMax',{'rescaleInv'});
            
        
   % CompareAllData(dataM,1);
   %CompareAllData(dataM,2);
   %CompareAllData(dataM,3);
    
    %dataM{l_diff}(Cak,y2Max)
    function PlotAsymptoticInterfaceResults(dataM,i_Cak,value,opts)
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
                Ca_k      = dataM{j0}(i_Cak,j2).config.optsPhys.Cak;     
                l_diff    = dataM{j0}(i_Cak,j2).config.optsPhys.l_diff;                
                col       = cols{mod(j0,nocols)+1};
                                
                if(IsOption(opts,'IsoInterface'))
                    val = dataM{j0}(i_Cak,j2).IsoInterface.(value); 
                    y  = dataM{j0}(i_Cak,j2).IsoInterface.y2/l_diff;
                elseif(IsOption(opts,'mu_y2'))
                    val = dataM{j0}(i_Cak,j2).mu_y2{value}.mu_0P*l_diff/Ca_k;
                    val = dataM{j0}(i_Cak,j2).mu_y2{value}.muP*l_diff^2/Ca_k;
                    y   = dataM{j0}(i_Cak,j2).mu_y2{value}.pts.y1_kv/l_diff;
                    y2  = dataM{j0}(i_Cak,j2).mu_y2{value}.pts.y2_kv(1)/l_diff;
                    %dataM(i,j).mu_y2{kk} = struct('mu_0',f_p,'pts',pts,'y2',y2(kk));
                end
                
                if(IsOption(opts,'SC'))
                    plot(y,val,['-',syms{mod(j2,nosyms)+1},col],...
                        'MarkerSize',8,'MarkerFaceColor',col);
                else
                    plot(y,val,['-',col],'linewidth',1.5); 
                end
                hold on;                                
                legendstr(end+1) = {['l_d = ',num2str(l_diff),' y_{2,Max}/l_d = ',num2str(dataM{j0}(i_Cak,j2).config.optsNum.PhysArea.y2Max/l_diff)]};
            end
        end  
        legend(legendstr);
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y_2/\ell_{d}$','Interpreter','Latex','fontsize',20);
        if(IsOption(opts,'mu_y2'))
            ylabel(['$\mu$ at $y_2/\ell_{d} = ',num2str(y2),'$'],'Interpreter','Latex','fontsize',20);  
        else
            ylabel(value,'fontsize',20);  
        end
    end        
    function PlotAsymptoticResults(dataM,parameter,opts)        
        
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
                    l_diff(j0)    = dataM{j0}(i_Cak,j2).config.optsPhys.l_diff;
                    if(IsOption(opts,'IsoInterface'))
                        par(j0)      = dataM{j0}(i_Cak,j2).IsoInterface.(parameter);                
                    else                        
                        par(j0)      = dataM{j0}(i_Cak,j2).(parameter);                
                    end
                end 
                col = cols{mod(j2,nocols)+1};
                if(IsOption(opts,'rescale'))
                    par = par./l_diff;
                elseif(IsOption(opts,'rescaleInv'))
                    par = par.*l_diff;
                end
                
                if(IsOption(opts,'rescaleCa'))
                    par = par/Ca;
                end
                plot(1./l_diff,par,...
                    ['-',syms{mod(i_Cak,nosyms)+1},col],...
                    'MarkerSize',8,'MarkerFaceColor',col); hold on;                
                
                legendstr(end+1) = {['Ca = ',num2str(Ca),' y_{2,Max}/l_d = ',num2str(dataM{1}(i_Cak,j2).config.optsNum.PhysArea.y2Max/l_diff(1))]};
             end  
        end
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        legend(legendstr);
        xlabel('$1/\ell_{d}$','Interpreter','Latex','fontsize',20);
        if(IsOption(opts,'rescale'))
            ylabel(['$',parameter,'/\ell_{d}$'],'Interpreter','Latex','fontsize',20);        
        else    
            ylabel(parameter,'Interpreter','Latex','fontsize',20);        
        end
        %ylim([90,100]);        
        SaveFigure([parameter,'_vs_l_diff']);
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
        for j = 1:length(pars.y2Max)            
            
            config.optsNum.PhysArea.y2Max = pars.y2Max(j);
            
            DI = DiffuseInterfaceBinaryFluid(config);
            DI.Preprocess();
            for i = 1:length(pars.Cak)
                
%                config.optsPhys.Cak = pars.Cak(i);
                
                DI.SetCak(pars.Cak(i));                
                
                DI.IterationStepFullProblem();                    
                
                DI.PostProcess(); 
                %DI.PlotResults();	  

                dataM(i,j).config      = DI.GetConfig();                
                dataM(i,j).Ca          = 3/4*pars.Cak(i);                                
                dataM(i,j).muMinusInf  = DI.mu(DI.IDC.Ind.left & DI.IDC.Ind.bottom);
                dataM(i,j).muPlusInf   = DI.mu(DI.IDC.Ind.right & DI.IDC.Ind.bottom);
                dataM(i,j).muMaxAbs    = max(abs(DI.mu));
                dataM(i,j).pMax        = max((DI.p));
                dataM(i,j).pMin        = min((DI.p));                              
                dataM(i,j).IsoInterface = DI.IsoInterface;                                
                dataM(i,j).stagnationPointY2 = DI.StagnationPoint.y2_kv(1);
                
                l_diff = DI.optsPhys.l_diff;
                
                y2 = 0:0.5:3;
                for kk = 1:length(y2)
                    [mu,pts]  = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.mu);
                    [muP,pts] = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy1*DI.mu);
                    dataM(i,j).mu_y2{kk} = struct('mu',mu,'muP',muP,'pts',pts,'y2',y2(kk));
                end
                
                close all;
            end
            clear('DI');
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
end