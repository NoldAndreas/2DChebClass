function ContactLineBinaryFluid
    close all;    
    AddPaths('DIBinaryPaper');            

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
    parameters.y2Max = 18:2:24;%(16:2:24);            
    parameters.l_d   = 1:0.25:3.0;%0.25:0.25:2.5;                    
	%parameters.l_d   = 2.25:0.25:3.0;%0.25:0.25:2.5;                    

    [dataM,~,res] = DataStorage('NumericalExperiment',@RunNumericalExperiment,parameters,[],[]);
    dataN = Rescale(dataM);clear('dataM');
    
    
    cols = {}; %{'g','b','c','k','r'};  
    nocols = 9;%length(cols);
    for iC = 1:nocols
        cols{end+1} = (nocols-iC)/nocols*[1 1 1];
	end
    %cols = {'g','b','c','k','r'};  nocols = length(cols);
	syms = {'<','d','s','o','>'};  nosyms = length(syms);            
    lines = {'-','--',':','.','-.'}; nolines = length(lines);                
    
    fileExampleVelocoties = 'D:\2DChebData\DIBinaryPaper\StagnationPoint_Velocity.fig';

    PlotAsymptoticInterfaceResults(dataN,1,[],struct('i',1,'val','flux1Slip'),{'mu_y2'});    
   % PlotExample(dataN,1,4,7,parameters);
                
    PlotAsymptoticResults(dataN,'hatL',{'IsoInterface','noLegend'}); xlim([0 1]); ylim([0.42 0.48]); SaveFigure('hatL_Asymptotics');
    
    
    f2 = open(fileExampleVelocoties); set(gca,'Xtick',[]);  set(gca,'Ytick',[]);
    PlotAsymptoticResults(dataN,'stagnationPointY2',{'noLegend'}); xlim([0 1]); ylim([2.4 2.55]);     
    inset2(gcf,f2,0.35,[0.25,0.65]); SaveFigure('stagnationPoint_Asymptotics');
    
    PlotAllInterfaceData(dataN,1,parameters);
    PlotAllInterfaceData(dataN,2,parameters);
   
    CompareAllData(dataN,1,parameters);
    CompareAllData(dataN,2,parameters);
   %CompareAllData(dataM,3,parameters);
    
    %Plot3D(dataN,1,4,parameters.config);
    PlotAsymptoticInterfaceResults(dataN,1,[],struct('i',1,'val','mu'),{'mu_y1'});    
    PlotAsymptoticInterfaceResults(dataN,1,[],struct('i',1,'val','mu_dy2'),{'mu_y2'});    
    
    PlotAllAsymptotics({'y2Max'},parameters);
    PlotAllAsymptotics({'Cn'},parameters);
        
    PlotAsymptoticResults(dataN,'d',{'mu_y2'});    
    SaveFigures(dataN,res,parameters);
     
  %  PlotAsymptoticInterfaceResults(dataN,1,'theta',{'IsoInterface'},'\theta');        
    PlotAsymptoticInterfaceResults(dataN,[],[],'kappa',{'IsoInterface'},'\kappa');    
    PlotAsymptoticInterfaceResults(dataN,1,[],'mu',{'IsoInterface'},'\mu');
    PlotAsymptoticInterfaceResults(dataN,1,[],'mu_ddy2',{'IsoInterface'},'\frac{\partial^2\mu}{\partial y_2^2}');
    PlotAsymptoticInterfaceResults(dataN,1,[],'p',{'IsoInterface'},'p');
    PlotAsymptoticInterfaceResults(dataN,1,[],'kappa_Cakmu',{'IsoInterface'},'\kappa/(Ca_k\mu)');        
	PlotAsymptoticInterfaceResults(dataN,1,[],'mu',{'mu_y2'},'\mu');        
    PlotAsymptoticInterfaceResults(dataN,1,[],'mu_ddy2',{'mu_y2'},'\frac{\partial^2\mu}{\partial y_2^2}');                                           
        
    function PlotExample(dataN,i_Cak,i_y2Max,i_Cn,parameters)
        data = dataN{i_Cn}(i_Cak,i_y2Max);
        
        config                        = parameters.config;
        Cn                            = data.Cn;  
        Cak                           = data.Cak;
        config.optsPhys.l_diff        = 1/Cn;
        config.optsPhys.Cak           = Cak;
        config.optsNum.PhysArea.y2Max = data.y2Max/Cn;   
        
        xL = [-7 7]; yL = [0 14];
        config.optsNum.PlotArea = struct('y1Min',xL(1)/Cn,'y1Max',xL(2)/Cn,...
                                         'y2Min',yL(1)/Cn,'y2Max',yL(2)/Cn,'N1',80,'N2',80);   

        DI = DiffuseInterfaceBinaryFluid(config);
        DI.Preprocess();
        DI.IterationStepFullProblem();                    
        DI.PostProcess();        
        
        PtsCart = DI.IDC.GetCartPts;
        
        %Plot Values at stagnation point
        y2_SP = DI.StagnationPoint.y2_kv(1);
        y1    = [-10,10]/Cn+DI.StagnationPoint.y1_kv(1);
        [mu,pts] = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.mu); close all; mu = mu/(Cn*Cak);
        [u,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.uv(1:end/2)); close all;
        [v,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.uv(1+end/2:end)); close all;
        [p,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.p); close all; p = p/Cn;
        [j1,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.flux(1:end/2)); close all;
        [j2,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.flux(1+end/2:end)); close all;
        
        lw  = 1.5;
        y   = (pts.y1_kv - DI.StagnationPoint.y1_kv(1))*Cn;
        figure('Position',[0 0 800 800],'color','white');
        plot(y,mu,'-k','linewidth',lw); hold on;
        plot(y,p,'--k','linewidth',lw);
        plot(y,u,'-b','linewidth',lw);
        plot(y,v,':b','linewidth',lw);
        plot(y,j1,'-m','linewidth',lw);
        plot(y,j2,':m','linewidth',lw);
        
        set(gca,'linewidth',1.5); set(gca,'fontsize',20);
        xlabel('$y_1$','fontsize',20,'Interpreter','Latex');       
        pbaspect([1 1 1]);        
        config.stagnationPoint = DI.StagnationPoint;
        SaveFigure('ValuesThroughStagnationPoint',config);               
                        
        %Plot velocities
        figure('Position',[0 0 800 800],'color','white');       
        DI.PlotResultsPhi();
        DI.PlotU();
        DI.PlotStagnationPoint();
       % DI.IDC.plotFlux(DI.uv,[],[],1,'k');%[y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] = 
       % DI.IDC.plotFlux(DI.flux); %[y1_s,y2_s,fl_y1_q,fl_y2_q] =         
        SetTicksLabels(Cn);                 
        SaveFigure('StagnationPoint_Velocity',config);
        
        %Plot flux j_1-j_2 next to stagnation point
        figure('Position',[0 0 800 800],'color','white');       
        DI.PlotResultsPhi();
        DI.PlotU(DI.flux); %[y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] = 
        DI.PlotStagnationPoint();
        %DI.IDC.plotFlux(DI.flux); %[y1_s,y2_s,fl_y1_q,fl_y2_q] =         
        SetTicksLabels(Cn);                 
        SaveFigure('StagnationPoint_Flux',config);
                
        %******* Plot 3D chemical potential, pressure and phi
        vals = {'mu','p','phi'};
        labs = {'$\mu$','p','$\phi$'};
        for k = 1:3            
            figure('Position',[0 0 800 800],'color','white');        
            if(strcmp(vals{k},'mu'))
                 DI.IDC.plot(DI.(vals{k})/(Cn*Cak)  );
            elseif(strcmp(vals{k},'p'))
                 DI.IDC.plot(DI.(vals{k})/Cn  );
            else
                 DI.IDC.plot(DI.(vals{k}));                 
            end                        
            DI.PlotU([],struct('color','b','linewidth',1.4));%close all;[y1M,y2M,VIM] =   %close all;[y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] =
            SetTicksLabels(Cn);
           % surf(y1M*Cn,y2M*Cn,VIM);  hold on;
           % h = streamline(y1MU*Cn,y2MU*Cn,fl_y1,fl_y2,startMask1*Cn,startMask2*Cn);   
           % xlim(xL); ylim(yL);
           % pbaspect([(xL(2)-xL(1))/(yL(2)-yL(1)) 1 1]);
           % set(h,'linewidth',1.5);  %set(h,'color','k');
            %set(gca,'linewidth',1.5); set(gca,'fontsize',20);
            %xlabel('$y_1$','fontsize',20,'Interpreter','Latex');
            %ylabel('$y_2$','fontsize',20,'Interpreter','Latex');
            zlabel(labs{k},'fontsize',25,'Interpreter','Latex');  
            SaveFigure([vals{k},'_3DPlot'],config);
        end  
        %         %Plot velocities next to stagnation point
%         [y1M,y2M,VIM] = DI.IDC.plot(DI.phi,'contour'); close all;
%         [y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] = DI.PlotU(); close all;                
%         [y1_s,y2_s,fl_y1_q,fl_y2_q] = DI.IDC.plotFlux(DI.uv); close all;
%         
% 
%         figure('Position',[0 0 800 800],'color','white');
%         [C,h] = contour(y1M*Cn,y2M*Cn,VIM,'linewidth',2.5);  hold on;        
%         h = streamline(y1MU*Cn,y2MU*Cn,fl_y1,fl_y2,startMask1*Cn,startMask2*Cn);   hold on;
%         quiver(y1_s*Cn,y2_s*Cn,fl_y1_q,fl_y2_q);        
%         %quiver(y1_s_Sepp,y2_s_Sepp,fl_y1_q__Sepp,fl_y2_q__Sepp,'color','m');        
%         plot(DI.IsoInterface.h*Cn,DI.IDC.Pts.y2*Cn,'k','linewidth',3);
%         plot(DI.StagnationPoint.y1_kv(1)*Cn,DI.StagnationPoint.y2_kv(1)*Cn,'o','MarkerFaceColor','m','MarkerSize',12)
%         
%         set(gca,'linewidth',1.5); set(gca,'fontsize',20);
%         xlabel('$y_1$','fontsize',20,'Interpreter','Latex');
%         ylabel('$y_2$','fontsize',20,'Interpreter','Latex');
%         pbaspect([1 1 1]);
%         SaveFigure('StagnationPoint',config);

        %************************************
        %************************************
        %************************************
        %Plot wider area
        DI.optsNum.PlotArea = struct('y1Min',-10/Cn,'y1Max',10/Cn,...
                                     'y2Min',0/Cn,'y2Max',20/Cn,'N1',80,'N2',80);   
        DI.IDC.InterpolationPlotCart(DI.optsNum.PlotArea,true);
        
         %Plot flux j_1-j_2 next to stagnation point
         %Plot velocities next to stagnation point         
         uSepp = GetSeppecherSolutionCart_Blurred([PtsCart.y1_kv - DI.deltaX,...
                                                PtsCart.y2_kv],1,0,0,DI.theta);
        [y1_s_Sepp,y2_s_Sepp,fl_y1_q__Sepp,fl_y2_q__Sepp] = DI.IDC.plotFlux(uSepp); close all;                                                             
        [y1MU_Sepp,y2MU_Sepp,fl_y1_Sepp,fl_y2_Sepp,startMask1_Sepp,startMask2_Sepp] = DI.PlotU(uSepp); close all;                
        
        [y1M,y2M,VIM] = DI.IDC.plot(DI.phi,'contour'); close all;
        [y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] = DI.PlotU(); close all;                
        [y1_s,y2_s,fl_y1_q,fl_y2_q] = DI.IDC.plotFlux(DI.uv); close all;
        

        figure('Position',[0 0 800 800],'color','white');
        [C,h] = contour(y1M*Cn,y2M*Cn,VIM,'linewidth',2.5);  hold on;        
        streamline(y1MU*Cn,y2MU*Cn,fl_y1,fl_y2,startMask1*Cn,startMask2*Cn);   hold on;
        h = streamline(y1MU_Sepp*Cn,y2MU_Sepp*Cn,fl_y1_Sepp,fl_y2_Sepp,startMask1_Sepp*Cn,startMask2_Sepp*Cn);   hold on;
        set(h,'color','m');
        quiver(y1_s*Cn,y2_s*Cn,fl_y1_q,fl_y2_q);        
        quiver(y1_s_Sepp*Cn,y2_s_Sepp*Cn,fl_y1_q__Sepp,fl_y2_q__Sepp,'color','m');        
        plot(DI.IsoInterface.h*Cn,DI.IDC.Pts.y2*Cn,'k','linewidth',3);
        plot(DI.StagnationPoint.y1_kv(1)*Cn,DI.StagnationPoint.y2_kv(1)*Cn,'o','MarkerFaceColor','m','MarkerSize',12)
        
        set(gca,'linewidth',1.5); set(gca,'fontsize',20);
        xlabel('$y_1$','fontsize',20,'Interpreter','Latex');
        ylabel('$y_2$','fontsize',20,'Interpreter','Latex');
        pbaspect([1 1 1]);
        SaveFigure('CoxComparison',config);                  
        
    end    
    function Plot3D(dataN,i_Cak,i_y2Max)        
        
        Ca = dataN{1}(i_Cak,i_y2Max).Ca;            
        title(['Ca = ',num2str(Ca)]);
        for j0 = 1:length(dataN)
            col  = cols{mod(j0-1,nocols)+1};         
            data = dataN{j0}(i_Cak,i_y2Max);
            surf(data.Pts.y1,data.Pts.y2,reshape(data.mu,data.Pts.N2,data.Pts.N1));
            %colormap(col);
            %Cn    = dataM{j0}(i_Cak,j2).Cn;             
            %DI.IDC.plot(dataN{j0}(i_Cak,i_y2Max).mu);  
            hold on;          
        end         
        xlim([-7.5 7.5]);
        ylim([0 15]);
    end    
    function PlotAllInterfaceData(dataN,i_Cak,parameters)
                       
        i_y2Max = 1:4;   
        xL      = [0 dataN{1}(1,i_y2Max(end)).y2Max];
        defaultOpts = {'IsoInterface','noLegend','noNewFigure'};
        figure('Position',[0 0 1000 1400],'color','white');
        
        subplot(4,2,1);        
        PlotAsymptoticInterfaceResults(dataN,i_Cak,i_y2Max,'kappa',defaultOpts);    
        xlim(xL);
        
        subplot(4,2,2);
        PlotAsymptoticInterfaceResults(dataN,i_Cak,i_y2Max,'mu',defaultOpts);
        xlim(xL);
        
        subplot(4,2,3);        
        %PlotAsymptoticInterfaceResults(dataN,1,'mu_ddy2',{'IsoInterface'},'\frac{\partial^2\mu}{\partial y_2^2}');
        PlotAsymptoticInterfaceResults(dataN,i_Cak,i_y2Max,'p',defaultOpts);        
        xlim(xL);
        
        subplot(4,2,4);
        PlotAsymptoticInterfaceResults(dataN,i_Cak,i_y2Max,'kappa_Cakmu',defaultOpts);     hold on;
        plot([0,24],[1.5 1.5],'k--','linewidth',1.5);
        ylim([1.2 1.8]);                
        xlim(xL);
        
        subplot(4,2,5);
        PlotAsymptoticInterfaceResults(dataN,i_Cak,i_y2Max,'u_n',defaultOpts);    
        plot([0,24],[0 0],'k--','linewidth',1.5);       
        ylim([-0.5 1.5]);
        xlim(xL);
        
        subplot(4,2,6);
        PlotAsymptoticInterfaceResults(dataN,i_Cak,i_y2Max,'u_t',defaultOpts);    
        plot([0,24],[0 0],'k--','linewidth',1.5);        
        ylim([-0.02 0.06]);                
        xlim(xL);
                
        subplot(4,2,7);
        PlotAsymptoticInterfaceResults(dataN,i_Cak,i_y2Max,'flux_n',defaultOpts);    
        plot([0,24],[0 0],'k--','linewidth',1.5);        
        xlim(xL);
        
        subplot(4,2,8);
        PlotAsymptoticInterfaceResults(dataN,i_Cak,i_y2Max,'flux_t',defaultOpts);    hold on;
        plot([0,24],[0 0],'k--','linewidth',1.5);
        xlim(xL);
        
        SaveFigure(['InterfaceData_Cak',num2str(dataN{1}(i_Cak,1).Cak)],parameters);                        
    end
    function PlotAllAsymptotics(opts,parameters)
        if(IsOption(opts,'y2Max'))
            plot_func = @PlotAsymptoticResults_Y2Max;
            name = 'y2Max';
        elseif(IsOption(opts,'Cn'))
            plot_func = @PlotAsymptoticResults;
            name = 'Cn';
        end
        
        defaultOpts = {'noLegend','noNewFigure'};
        
        figure('Position',[0 0 1500 1500],'color','white');
        subplot(3,3,1);
        plot_func(dataN,'muMinusInf',defaultOpts);
        subplot(3,3,2);        
        plot_func(dataN,'muPlusInf',defaultOpts);
        subplot(3,3,3);
        plot_func(dataN,'stagnationPointY2',defaultOpts)
        subplot(3,3,4);
        plot_func(dataN,'hatL',[{'IsoInterface'},defaultOpts]);
        subplot(3,3,5);
        plot_func(dataN,'muMaxAbs',defaultOpts);
        subplot(3,3,6);
        plot_func(dataN,'pMin',defaultOpts);
        subplot(3,3,7);
        plot_func(dataN,'pMax',defaultOpts);
        subplot(3,3,8);
        plot_func(dataN,'d',[{'mu_y2'},defaultOpts]);
        
        SaveFigure(['AsymptoticResults_',name],parameters);
    end
    function SaveFigures(dataN,res,params) 
        params.comment = 'run by ContactLineBinaryFluid.m';
        [~,fn]   = fileparts(res.Filename);        
        filename = ['NumericalExperiment' filesep fn '_'];
        
        PlotAsymptoticResults(dataN,'mu_max_y20',{'mu_y2','noLegend'});    
        SaveFigure([filename, 'wall_max_mu'],params);      
        
        PlotAsymptoticResults(dataN,'mu_dy1_max_y20_sqrtCn',{'mu_y2','noLegend'});    
        SaveFigure([filename, 'wall_max_mudy1_sqrtCn'],params);        
        
        PlotAsymptoticResults(dataN,'mu_ddy2_max_y20_sqrtCn',{'mu_y2','noLegend'});
        SaveFigure([filename, 'wall_max_muddy2_sqrtCn'],params);        
                
        PlotAsymptoticInterfaceResults(dataN,[],[],'mu',{'mu_y2','noLegend'},'\mu');                
        SaveFigure([filename, 'wall_mu'],params);
        
        PlotAsymptoticInterfaceResults(dataN,[],[],'mu_dy1',{'mu_y2','noLegend'},'\frac{\partial\mu}{\partial y_1}');                        
        SaveFigure([filename, 'wall_mu_dy1'],params);
        
        PlotAsymptoticInterfaceResults(dataN,[],[],'mu_ddy2',{'mu_y2','noLegend'},'\frac{\partial^2\mu}{\partial y_2^2}');                        
        SaveFigure([filename, 'wall_mu_ddy2'],params);
        
	    
        PlotAsymptoticResults(dataN,'hatL',{'IsoInterface','noLegend'}); 
        SaveFigure([filename, 'hatL'],params);
                       
        PlotAsymptoticResults(dataN,'stagnationPointY2',{'noLegend'});         
        SaveFigure([filename, 'stagnationPointY2'],params);
              
        SaveFigure([filename, 'Interface_KappaCak_Mu'],params);
    end
    function PlotAsymptoticInterfaceResults(dataM,i_Cak,i_y2Max,value,opts)
        if(nargin < 4)
            opts = {};                        
        end        
        if(ischar(opts))
            opts = {opts};
        end
        legendstr = {};
        
        if(isempty(i_Cak))
            i_Cak =  1:size(dataM{1},1);
        end
        if(isempty(i_y2Max))
            i_y2Max = size(dataM{1},2);
        end
        
        ylabelStr = GetYStr(value,opts,dataM);
            
        if(~IsOption(opts,'noNewFigure'))
            figure('Position',[0 0 800 800],'color','white');        
        end
        if(length(i_Cak) == 1)
            title(['Ca = ',num2str(dataM{1}(i_Cak,1).Ca)]);
        end
        
        for j1 = i_Cak
            lin   = lines{mod(j1-1,nolines)+1};
            for j2 = i_y2Max  
                sym = syms{mod(j2-1,nosyms)+1};                
                for j0 = 1:length(dataM)                
                    Cn    = dataM{j0}(j1,j2).Cn;                
                    col   = cols{mod(j0-1,nocols)+1};

                    if(IsOption(opts,'IsoInterface'))
                        val = dataM{j0}(j1,j2).IsoInterface.(value); 
                        y  = dataM{j0}(j1,j2).IsoInterface.y2;
                    elseif(IsOption(opts,'mu_y2'))                        
                        val = dataM{j0}(j1,j2).mu_y2{value.i}.(value.val);                        
                        y  = dataM{j0}(j1,j2).mu_y2{value.i}.pts.y1_kv - dataN{j0}(j1,j2).mu_y2{value.i}.y1_contactLine;   
                    elseif(IsOption(opts,'mu_y1'))
                        val = dataM{j0}(j1,j2).mu_y1{value.i}.(value.val);
                        y  = dataM{j0}(j1,j2).mu_y1{value.i}.pts.y2_kv;   
                    end

                    if(IsOption(opts,'SC'))
                        plot(y,val,[lin,sym],...
                            'color',col,...
                            'MarkerSize',8,'MarkerFaceColor',col);
                    else
                        plot(y,val,[lin],'color',col,'linewidth',1.5); 
                    end
                    hold on;                                
                    if(size(dataM{1},2)>1)
                        legendstr(end+1) = {['Cn = ',num2str(Cn),' y_{2,Max} = ',num2str(dataM{j0}(j1,j2).y2Max)]};
                    else
                        legendstr(end+1) = {['Cn = ',num2str(Cn)]};
                    end
                end
            end  
        end
        if(~IsOption(opts,'noLegend'))            
            legend(legendstr);
        end
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);        
        
        if(IsOption(opts,'mu_y2'))
            xlabel('$y_1 - y_{1,CL}$','Interpreter','Latex','fontsize',20);
        else
            xlabel('$y_2$','Interpreter','Latex','fontsize',20);
        end
        ylabel(ylabelStr,'Interpreter','Latex','fontsize',20);  
                
    end        
    function PlotAsymptoticResults(dataM,parameter,opts)        
        
        if(nargin == 2)
            opts = {};                        
        end        
        if(ischar(opts))
            opts = {opts};
        end                
        ylabelStr = GetYStr(parameter);
        
        
        if(~IsOption(opts,'noNewFigure'))            
            figure('Position',[0 0 800 800],'color','white');
        end
        legendstr = {}; 
        for i_Cak = 1:size(dataM{1},1)
            Ca = dataM{1}(i_Cak,1).Ca;
            for j2 = 1:size(dataM{1},2)
                for j0 = 1:length(dataM)
                    Cn(j0)    = dataM{j0}(i_Cak,j2).Cn;
                    if(IsOption(opts,'IsoInterface'))
                        par(j0)      = dataM{j0}(i_Cak,j2).IsoInterface.(parameter);                
                    elseif(IsOption(opts,'mu_y2'))
                        par(j0)      = dataM{j0}(i_Cak,j2).mu_y2{1}.(parameter);                                        
                    else                        
                        par(j0)      = dataM{j0}(i_Cak,j2).(parameter);                
                    end
                end                 
                %col  = cols{mod(j0-1,nocols)+1};                                
                lin = lines{mod(i_Cak-1,nolines)+1};                
                sym = syms{mod(j2-1,nosyms)+1};
                                
                plot(Cn,par,...
                    [lin,sym,'k'],'linewidth',1.5,...
                    'MarkerSize',8,'MarkerFaceColor','k'); hold on;                
                
                legendstr(end+1) = {['Ca = ',num2str(Ca),' y_{2,Max} = ',num2str(dataM{1}(i_Cak,j2).y2Max)]};
             end  
        end        
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',25);
        if(~IsOption(opts,'noLegend'))            
            legend(legendstr);
        end
        xlabel('Cn','Interpreter','Latex','fontsize',25);        
        ylabel(ylabelStr,'Interpreter','Latex','fontsize',25);                
        
        %ylim([90,100]);        
        %SaveFigure([parameter,'_vs_l_diff']);
    end        
    function PlotAsymptoticResults_Y2Max(dataM,parameter,opts)
        if(nargin == 2)
            opts = {};                        
        end        
        if(ischar(opts))
            opts = {opts};
        end                
        
        if(~IsOption(opts,'noNewFigure'))            
            figure('Position',[0 0 800 800],'color','white');
        end
         
        legendstr = {};
        for j1 = 1:size(dataM{1},1)
            Ca = dataM{1}(j1,1).Ca;
            for j0 = 1:length(dataM)                
                Cn = dataM{j0}(1,1).Cn;
                for j2 = 1:size(dataM{1},2)
                    y2M(j2)  = dataM{j0}(j1,j2).y2Max;                                        
                    if(IsOption(opts,'IsoInterface'))
                        par(j2)      = dataM{j0}(j1,j2).IsoInterface.(parameter);                
                    elseif(IsOption(opts,'mu_y2'))
                        par(j2)      = dataM{j0}(j1,j2).mu_y2{1}.(parameter);                                        
                    else                        
                        par(j2)      = dataM{j0}(j1,j2).(parameter);                
                    end                                        
                end                            
                col  = cols{mod(j0-1,nocols)+1};                
                lin = lines{mod(j1-1,nolines)+1};                
                %sym = syms{mod(j2,nosyms)+1};
                
                plot(y2M,par,...
                    lin,'color',col,'linewidth',1.5,...
                    'MarkerSize',8,'MarkerFaceColor',col); hold on;                
                
                legendstr(end+1) = {['Ca = ',num2str(Ca),' Cn = ',num2str(Cn)]};
             end  
        end
        
        if(~IsOption(opts,'noLegend'))            
            legend(legendstr);
        end
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y_{2,max}$','Interpreter','Latex','fontsize',20);       
        ylabel(GetYStr(parameter),'Interpreter','Latex','fontsize',20);        
        %ylim([90,100]);
        %SaveFigure([parameter,'_vs_y2Max']);
    end
    function CompareAllData(dataM,j1,parameters)
        thetaEq = dataM{1}(1,1).thetaEq;        

        figure('Position',[0 0 800 800],'color','white');
        hatL_M = zeros(size(dataM));
        Sy2_M  = zeros(size(dataM));
        Ca     = dataM{1}(j1,1).Ca;
        
        for j0 = 1:length(dataM)
            for j2 = 1:size(dataM{j0},2)
                hatL_M(j1,j2) = dataM{j0}(j1,j2).IsoInterface.hatL;         
                Sy2_M(j1,j2)  = dataM{j0}(j1,j2).stagnationPointY2; 
                
                col  = cols{mod(j0-1,nocols)+1};                                
                %lin = lines{mod(i_Cak-1,nolines)+1};                
                sym = syms{mod(j2-1,nosyms)+1};
                
                PlotInterfaceSlope(dataM{j0},j1,j2,col,sym);
            end                                       
        end       
                
        y2P       = (1:0.1:max(dataM{end}(j1,j2).IsoInterface.y2 - 4))';
        hatL      = dataM{end}(j1,j2).IsoInterface.hatL;
        theta_Ana = GHR_Inv(Ca*log(y2P/hatL)+GHR_lambdaEta(thetaEq,1),1);
        plot(y2P,180/pi*theta_Ana,'m','linewidth',3); hold on;        


        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y$','Interpreter','Latex','fontsize',20);
        ylabel('$\theta[^\circ]$','Interpreter','Latex','fontsize',20);        
        %ylim([90,100]);
                
        SaveFigure(['InterfaceSlope_Ca_',num2str(Ca)],parameters);
        
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
    function PlotInterfaceSlope(dataM,i1,i2,col,sym)                    
        y2     = dataM(i1,i2).IsoInterface.y2;                                          
        mark   = (y2 < (max(y2)-4));
        plot(y2(mark),180/pi*dataM(i1,i2).IsoInterface.theta(mark),['-',sym],...
            'color',col,'MarkerSize',7,'MarkerFaceColor',col); hold on;                                              
    end
    function dataN = Rescale(dataM)
        for k0 = 1:length(dataM)
            for k1 = 1:size(dataM{1},1)
                for k2 = 1:size(dataM{1},2)
                    Cn    = 1/dataM{k0}(k1,k2).config.optsPhys.l_diff;
                    Cak   = dataM{k0}(k1,k2).config.optsPhys.Cak;

                    dataN{k0}(k1,k2).Cn          = Cn;
                    dataN{k0}(k1,k2).Cak         = Cak;
                    dataN{k0}(k1,k2).y2Max       = dataM{k0}(k1,k2).config.optsNum.PhysArea.y2Max*Cn;
                    dataN{k0}(k1,k2).Ca          = dataM{k0}(k1,k2).Ca;
                    dataN{k0}(k1,k2).thetaEq     = dataM{k0}(k1,k2).config.optsPhys.thetaEq;

                    dataN{k0}(k1,k2).Pts.N1      = dataM{k0}(k1,k2).Pts.N1;
                    dataN{k0}(k1,k2).Pts.N2      = dataM{k0}(k1,k2).Pts.N2;
                    dataN{k0}(k1,k2).Pts.y1      = dataM{k0}(k1,k2).Pts.y1*Cn;                                
                    dataN{k0}(k1,k2).Pts.y2      = dataM{k0}(k1,k2).Pts.y2*Cn;                                
                    dataN{k0}(k1,k2).Pts.y1_kv   = dataM{k0}(k1,k2).Pts.y1_kv*Cn;                                
                    dataN{k0}(k1,k2).Pts.y2_kv   = dataM{k0}(k1,k2).Pts.y2_kv*Cn; 
                    
                    dataN{k0}(k1,k2).muMinusInf  = dataM{k0}(k1,k2).muMinusInf/(Cn*Cak);
                    dataN{k0}(k1,k2).muPlusInf   = dataM{k0}(k1,k2).muPlusInf/(Cn*Cak);
                    dataN{k0}(k1,k2).muMaxAbs    = dataM{k0}(k1,k2).muMaxAbs/(Cn*Cak);
                    dataN{k0}(k1,k2).pMax        = dataM{k0}(k1,k2).pMax/Cn;
                    dataN{k0}(k1,k2).pMin        = dataM{k0}(k1,k2).pMin/Cn;                                      
                    dataN{k0}(k1,k2).stagnationPointY2 = dataM{k0}(k1,k2).stagnationPointY2*Cn;
                    dataN{k0}(k1,k2).mu          = dataM{k0}(k1,k2).mu/(Cn*Cak);

                    for i = 1:2
                        if(i==1)
                            mu_str = 'mu_y2';
                        elseif(i==2)
                            mu_str = 'mu_y1';
                        end
                            
                        for kk = 1:length(dataM{1}(1,1).(mu_str))
                            dataN{k0}(k1,k2).(mu_str){kk}.mu        = dataM{k0}(k1,k2).(mu_str){kk}.mu/(Cn*Cak);
                            dataN{k0}(k1,k2).(mu_str){kk}.mu_dy1    = dataM{k0}(k1,k2).(mu_str){kk}.mu_dy1/(Cn^2*Cak);
                            dataN{k0}(k1,k2).(mu_str){kk}.mu_dy2    = dataM{k0}(k1,k2).(mu_str){kk}.mu_dy2/(Cn^2*Cak);
                            dataN{k0}(k1,k2).(mu_str){kk}.mu_ddy2   = dataM{k0}(k1,k2).(mu_str){kk}.mu_ddy2/(Cn^3*Cak);
                            dataN{k0}(k1,k2).(mu_str){kk}.pts.y1_kv = dataM{k0}(k1,k2).(mu_str){kk}.pts.y1_kv*Cn;
                            dataN{k0}(k1,k2).(mu_str){kk}.pts.y2_kv = dataM{k0}(k1,k2).(mu_str){kk}.pts.y2_kv*Cn;
                            
                            if(strcmp(mu_str,'mu_y2'))
                                dataN{k0}(k1,k2).(mu_str){kk}.y2        = dataM{k0}(k1,k2).(mu_str){kk}.y2*Cn;
                                dataN{k0}(k1,k2).(mu_str){kk}.flux1     = dataM{k0}(k1,k2).(mu_str){kk}.flux1;
                                dataN{k0}(k1,k2).(mu_str){kk}.flux1Slip = 1 - abs(dataM{k0}(k1,k2).(mu_str){kk}.flux1);
                                dataN{k0}(k1,k2).(mu_str){kk}.phi       = dataM{k0}(k1,k2).(mu_str){kk}.phi;
                                dataN{k0}(k1,k2).(mu_str){kk}.phiflux1  = (dataM{k0}(k1,k2).(mu_str){kk}.flux1)./(dataM{k0}(k1,k2).(mu_str){kk}.phi);
                                
                                %Get Contact Line Position:
                                y1  = dataN{k0}(k1,k2).(mu_str){kk}.pts.y1_kv;
                                phi = dataN{k0}(k1,k2).(mu_str){kk}.phi;
                                i1 = find(phi < 0,1,'last');
                                i2 = find(phi > 0,1,'first');
                                if((i2-i1) ~= 1)
                                    error('something weird happened');
                                end
                                y10 = y1(i1) - (y1(i2)-y1(i1))*phi(i1)/(phi(i2)-phi(i1));
                                dataN{k0}(k1,k2).(mu_str){kk}.y1_contactLine = y10;

                            elseif(strcmp(mu_str,'mu_y1'))
                                dataN{k0}(k1,k2).(mu_str){kk}.y1        = dataM{k0}(k1,k2).(mu_str){kk}.y1*Cn;
                            end

                            dataN{k0}(k1,k2).(mu_str){kk}.mu_max_y20             = max(abs(dataN{k0}(k1,k2).(mu_str){kk}.mu));
                            [muMax,i] = max((dataN{k0}(k1,k2).(mu_str){kk}.mu_dy1));
                            [muMin,j] = min((dataN{k0}(k1,k2).(mu_str){kk}.mu_dy1));
                            dataN{k0}(k1,k2).(mu_str){kk}.d                      = dataN{k0}(k1,k2).(mu_str){kk}.pts.y1_kv(i)-dataN{k0}(k1,k2).(mu_str){kk}.pts.y1_kv(j);
                            dataN{k0}(k1,k2).(mu_str){kk}.mu_dy1_max_y20_sqrtCn  = max(abs(dataN{k0}(k1,k2).(mu_str){kk}.mu_dy1))*sqrt(Cn);
                            dataN{k0}(k1,k2).(mu_str){kk}.mu_ddy2_max_y20_sqrtCn = max(abs(dataN{k0}(k1,k2).(mu_str){kk}.mu_ddy2))*sqrt(Cn);
                        end
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
    function str = GetYStr(parameter,opts,data)
        
        if((nargin > 1) && (IsOption(opts,'mu_y1') || IsOption(opts,'mu_y2')))
            VarName = parameter.val;
        else
            VarName = parameter;
        end
        
        if(strcmp(VarName,'d'))
            str = '$d$ at $y_2 = 0$';        
        elseif(strcmp(VarName,'flux_t'))
            str = '$j_t$';
        elseif(strcmp(VarName,'flux_n'))
            str = '$j_n$';
        elseif(strcmp(VarName,'u_t'))
            str = '$u_t$';
        elseif(strcmp(VarName,'u_n'))
            str = '$u_n$';
        elseif(strcmp(VarName,'kappa_Cakmu'))
            str = '$\kappa/(Ca_k\mu)$';
        elseif(strcmp(VarName,'kappa'))
            str = '$\kappa$';
        elseif(strcmp(VarName,'mu'))
            str = '$\mu$';
        elseif(strcmp(VarName,'mu_dy1'))
            str = '$\frac{\partial \mu}{\partial y_1}$';
        elseif(strcmp(VarName,'mu_dy2'))
            str = '$\frac{\partial \mu}{\partial y_2}$';
        elseif(strcmp(VarName,'muMinusInf'))
            str = '$\mu_{y_1 = -\infty}$';
        elseif(strcmp(VarName,'muPlusInf'))
            str = '$\mu_{y_1 = \infty}$';
        elseif(strcmp(VarName,'stagnationPointY2'))
            str = '$y_{2,S}$';
        elseif(strcmp(VarName,'hatL'))
            str = '$\hat L$';
        elseif(strcmp(VarName,'muMaxAbs'))
            str = '$\max|\mu|$';   
        elseif(strcmp(VarName,'pMin'))
            str = '$\min p$';            
        elseif(strcmp(VarName,'pMax'))
            str = '$\max p$';           
        elseif(strcmp(VarName,'mu_max_y20'))
            str = '$\max|\mu|$ at $y_2 = 0$';
        elseif(strcmp(VarName,'mu_dy1_max_y20_sqrtCn'))
            str = '$\sqrt{Cn}\max|\frac{\partial \mu}{\partial y_1}|$ at $y_2 = 0$';
        elseif(strcmp(VarName,'mu_ddy2_max_y20_sqrtCn'))
            str = '$\sqrt{Cn}\max|\frac{\partial^2 \mu}{\partial y_2^2}|$ at $y_2 = 0$';
        elseif(strcmp(VarName,'d'))
            str = '$d$ at $y_2 = 0$';
        else
            str = VarName;
        end
        
        if(nargin > 1)
            if(IsOption(opts,'mu_y1'))
                str = [str,' at $y_1 = ',num2str(data{1}(1,1).mu_y1{parameter.i}.y1),'$'];
            elseif(IsOption(opts,'mu_y2'))
                str = [str,' at $y_2 = ',num2str(data{1}(1,1).mu_y2{parameter.i}.y2),'$'];                          
            end
        end
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
                    
                    dataM{k0}(i,j).Pts = DI.IDC.Pts;
                    dataM{k0}(i,j).Diff.Dy1  = DI.IDC.Diff.Dy1;
                    dataM{k0}(i,j).Diff.Dy2  = DI.IDC.Diff.Dy2;
                    dataM{k0}(i,j).Diff.DDy2 = DI.IDC.Diff.DDy2;
                    dataM{k0}(i,j).Diff.DDy1 = DI.IDC.Diff.DDy1;
                    dataM{k0}(i,j).mu  = DI.mu;
                    dataM{k0}(i,j).uv  = DI.uv;
                    dataM{k0}(i,j).p   = DI.p;
                    dataM{k0}(i,j).phi = DI.phi;
                    
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
                        interval = l_diff*[-50,50];
                        [mu,pts]      = DI.IDC.plotLine(interval,l_diff*y2(kk)*[1 1],DI.mu);                    
                        [mu_dy1,pts]  = DI.IDC.plotLine(interval,l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy1*DI.mu);
                        [mu_dy2,pts]  = DI.IDC.plotLine(interval,l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy2*DI.mu);
                        [mu_ddy2,pts] = DI.IDC.plotLine(interval,l_diff*y2(kk)*[1 1],DI.IDC.Diff.DDy2*DI.mu);
                        [flux1,pts]   = DI.IDC.plotLine(interval,l_diff*y2(kk)*[1 1],DI.flux(1:end/2));                    
                        [phi,pts]     = DI.IDC.plotLine(interval,l_diff*y2(kk)*[1 1],DI.phi);            
                        %[v_dy2,pts] = DI.IDC.plotLine(interval,l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy2*DI.uv(1+end/2:end));
                        dataM{k0}(i,j).mu_y2{kk} = struct('mu',mu,...%'v_dy2',v_dy2,...
                                                          'phi',phi,...
                                                          'flux1',flux1,...
                                                          'mu_dy1',mu_dy1,...
                                                          'mu_dy2',mu_dy2,...
                                                          'mu_ddy2',mu_ddy2,...
                                                          'pts',pts,...
                                                          'y2',y2(kk));                        
                    end
                    
                    
                    y1d = (-2:0.5:3);                    
                    for kk = 1:length(y1d)
                        y1            = DI.StagnationPoint.y1_kv(1)+y1d(kk);
                        [mu,pts]      = DI.IDC.plotLine(l_diff*y1*[1 1],l_diff*[0 10],DI.mu);                    
                        [mu_dy1,pts]  = DI.IDC.plotLine(l_diff*y1*[1 1],l_diff*[0 10],DI.IDC.Diff.Dy1*DI.mu);
                        [mu_dy2,pts]  = DI.IDC.plotLine(l_diff*y1*[1 1],l_diff*[0 10],DI.IDC.Diff.Dy2*DI.mu);
                        [mu_ddy2,pts] = DI.IDC.plotLine(l_diff*y1*[1 1],l_diff*[0 10],DI.IDC.Diff.DDy2*DI.mu);
                        %[v_dy2,pts] = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy2*DI.uv(1+end/2:end));
                        dataM{k0}(i,j).mu_y1{kk} = struct('mu',mu,...%'v_dy2',v_dy2,...
                                                      'mu_dy1',mu_dy1,...
                                                      'mu_dy2',mu_dy2,...
                                                      'mu_ddy2',mu_ddy2,...
                                                      'pts',pts,'y1',y1,'y1d',y1d(kk));                        
                    end

                    close all;
                end
                clear('DI');
            end
        end        
    end  

    function SetTicksLabels(Cn)
        T = get(gca,'Xtick');                
        c = {};
        for i = 1:length(T)
            c{end+1} = num2str(T(i)*Cn);
        end
        set(gca,'XTickLabel',c);
        
        
        T = get(gca,'Ytick');                
        c = {};
        for i = 1:length(T)
            c{end+1} = num2str(T(i)*Cn);
        end
        set(gca,'YTickLabel',c);                          
        
        set(gca,'linewidth',1.5); set(gca,'fontsize',25);
        xlabel('$y_1$','fontsize',25,'Interpreter','Latex');
        ylabel('$y_2$','fontsize',25,'Interpreter','Latex');        
        pbaspect([1 1 1]);
     
    end
end


%     function dataM = RunNumericalExperiment(pars,h)
%     
%         config = pars.config;
%         for k0 = 1:length(pars.l_d)
%             
%             l_diff                 = pars.l_d(k0);
%             config.optsPhys.l_diff = l_diff;    
%             
%             for j = 1:length(pars.y2Max)            
% 
%                 config.optsNum.PhysArea.y2Max = pars.y2Max(j)*l_diff;
% 
%                 DI = DiffuseInterfaceBinaryFluid(config);
%                 DI.Preprocess();
%                 for i = 1:length(pars.Cak)
% 
%                     DI.SetCak(pars.Cak(i));                
%                     DI.IterationStepFullProblem();                    
%                     DI.PostProcess();                       
% 
%                     dataM{k0}(i,j).config.optsNum    = DI.optsNum;
%                     dataM{k0}(i,j).config.optsPhys   = DI.optsPhys;
%                     
%                     dataM{k0}(i,j).Pts = DI.IDC.Pts;
%                     dataM{k0}(i,j).mu  = DI.mu;
%                     dataM{k0}(i,j).uv  = DI.uv;
%                     dataM{k0}(i,j).p   = DI.p;
%                     dataM{k0}(i,j).phi = DI.phi;
%                     
%                     dataM{k0}(i,j).Ca                = 3/4*pars.Cak(i);                                
%                     dataM{k0}(i,j).muMinusInf        = DI.mu(DI.IDC.Ind.left & DI.IDC.Ind.bottom);
%                     dataM{k0}(i,j).muPlusInf         = DI.mu(DI.IDC.Ind.right & DI.IDC.Ind.bottom);
%                     dataM{k0}(i,j).muMaxAbs          = max(abs(DI.mu));
%                     dataM{k0}(i,j).pMax              = max((DI.p));
%                     dataM{k0}(i,j).pMin              = min((DI.p));                              
%                     dataM{k0}(i,j).IsoInterface      = DI.IsoInterface;                                
%                     dataM{k0}(i,j).stagnationPointY2 = DI.StagnationPoint.y2_kv(1);                    
% 
%                     y2 = 0:0.5:3;
%                     for kk = 1:length(y2)
%                         [mu,pts]      = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.mu);                    
%                         [mu_dy1,pts]  = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy1*DI.mu);
%                         [mu_dy2,pts]  = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy2*DI.mu);
%                         [mu_ddy2,pts] = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.DDy2*DI.mu);
%                         %[v_dy2,pts] = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy2*DI.uv(1+end/2:end));
%                         dataM{k0}(i,j).mu_y2{kk} = struct('mu',mu,...%'v_dy2',v_dy2,...
%                                                       'mu_dy1',mu_dy1,...
%                                                       'mu_dy2',mu_dy2,...
%                                                       'mu_ddy2',mu_ddy2,...
%                                                       'pts',pts,'y2',y2(kk));                        
%                     end
%                     
%                     
%                     y1 = -2:0.5:3;
%                     for kk = 1:length(y2)
%                         [mu,pts]      = DI.IDC.plotLine(l_diff*y1(kk)*[1 1],l_diff*[0 10],DI.mu);                    
%                         [mu_dy1,pts]  = DI.IDC.plotLine(l_diff*y1(kk)*[1 1],l_diff*[0 10],DI.IDC.Diff.Dy1*DI.mu);
%                         [mu_dy2,pts]  = DI.IDC.plotLine(l_diff*y1(kk)*[1 1],l_diff*[0 10],DI.IDC.Diff.Dy2*DI.mu);
%                         [mu_ddy2,pts] = DI.IDC.plotLine(l_diff*y1(kk)*[1 1],l_diff*[0 10],DI.IDC.Diff.DDy2*DI.mu);
%                         %[v_dy2,pts] = DI.IDC.plotLine(l_diff*[-10 10],l_diff*y2(kk)*[1 1],DI.IDC.Diff.Dy2*DI.uv(1+end/2:end));
%                         dataM{k0}(i,j).mu_y1{kk} = struct('mu',mu,...%'v_dy2',v_dy2,...
%                                                       'mu_dy1',mu_dy1,...
%                                                       'mu_dy2',mu_dy2,...
%                                                       'mu_ddy2',mu_ddy2,...
%                                                       'pts',pts,'y1',y1(kk));                        
%                     end
% 
%                     close all;
%                 end
%                 clear('DI');
%             end
%         end        
%     end  