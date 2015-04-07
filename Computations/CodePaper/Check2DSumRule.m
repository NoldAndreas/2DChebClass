function Check2DSumRule()
    
    AddPaths('CodePaper');            
    
    PhysArea = struct('N',[20,20],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',90);
                      
    V2Num    = struct('Fex','SplitDisk','N',[80,80]);
    V2       = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',5);     

    FexNum   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);

    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',FexNum,'V2Num',V2Num,...
                     'maxComp_y2',25,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.0);%1.375);%1.25)s;
    

    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);

    config = v2struct(optsNum,optsPhys);                        

    Interval_SumRule = [-10,20];
    N    = 20:10:60;        
    
    ignoreList = {'config_optsNum_PhysArea_N',...
                  'config_optsNum_PhysArea_N2bound'};
    comp = [];
    dirRes = 'SumRuleError2D';
    
    AddPaths('CodePaper');    
    config.optsNum.V2Num.Fex     = 'SplitAnnulus'; 
    config.optsPhys.V2.V2DV2     = 'BarkerHendersonCutoff_2D';                
    config.optsPhys.V1.epsilon_w = 0.94;                
    res{1}                       = DataStorage(dirRes,@ComputeError,v2struct(N,config,Interval_SumRule),[],comp,ignoreList);     
                
    config.optsNum.V2Num.Fex     = 'SplitDisk'; 
    config.optsPhys.V2.V2DV2     = 'ExponentialDouble';                          
    config.optsPhys.V1.V1DV1     = 'Vext_Exp_HardWall';
    config.optsPhys.V1.epsilon_w = 1.45;         
    res{2} = DataStorage(dirRes,@ComputeError,v2struct(N,config,Interval_SumRule),[],comp,ignoreList);                
        
    config.optsNum.V2Num.Fex     = 'SplitDisk';     
    config.optsPhys.V2.V2DV2     = 'BarkerHenderson_2D';     
    config.optsPhys.V1.V1DV1     = 'Vext_BarkerHenderson_HardWall';    
    config.optsPhys.V1.epsilon_w = 1.0;                
    res{3} = DataStorage(dirRes,@ComputeError,v2struct(N,config,Interval_SumRule),[],comp,ignoreList);
    
     %AddPaths('CodePaper');
     %config.optsNum.V2Num.Fex = 'SplitAnnulus';
     %config.optsPhys.V2.V2DV2 = 'BarkerHendersonHardCutoff_2D';                
%     res{2}                   = DataStorage(dirRes,@ComputeError,v2struct(N,config),[],comp,ignoreList);                        
    
    cols = GetGreyShades(length(res{1}));
	syms = {'o','^','*','<','d','s','>'};  nosyms = length(syms);
    lines = {'-','--',':','.','-.'}; nolines = length(lines);                

    
    %**********************
    %**********************
    PlotDisjoiningPressureProfiles(res);    
    %**********************
    %**********************
    %PlotDensityProfile(res{1}(3));
    %PlotDensityProfile(res{2}(3));
    %PlotDensityProfile(res{3}(3));         
    %**********************
    %**********************
    figure('color','white','Position',[0 0 2000 800]); 
    legendstring = {};
    subplot(1,2,1);
    PlotErrorGraph(res,'sumRuleII_relError','k');            
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);            
    xlabel('$N$','Interpreter','Latex','fontsize',20);
    ylabel(['Relative sum rule error $\frac{\int_{-\infty}^\infty\Pi_{II} dx + \gamma_{lv}\sin \theta_{{Y}}}{\gamma_{lv}\sin \theta_Y}$'],...
                            'Interpreter','Latex','fontsize',20);
    xlim([(N(1)-2),(N(end)+2)]);    
    
    subplot(1,2,2);
    PlotErrorGraph(res,'sumRuleII_eps','k','eps');     
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);            
    xlabel('$N$','Interpreter','Latex','fontsize',20);
    ylabel('$\frac{\Delta \Pi_{II}}{\max|\Pi_{II}|}$',...
                            'Interpreter','Latex','fontsize',20);
    xlim([(N(1)-2),(N(end)+2)]);    
    %legend(legendstring,'Location','northeast');%,'Orientation','horizontal');
   % legend(legendstring,'Location','eastOutside');%,'Orientation','horizontal');
        
    comment = '';
    SaveFigure('SumRuleError2D',v2struct(N,config,comment));
    
    function PlotDisjoiningPressureProfiles(res)
        
        K = length(res);
        y = res{1}(1).piII.y1.Pts.y;
        figure('color','white','Position',[0 0 2000 500]); 
        
        for k0 = 1:K
            subplot(1,K,k0);
            resK = res{k0};
            for i = 1:length(resK)
                plot(y,...
                     resK(i).piII.DP,['-',syms{k0}],'color',cols{i},...
                                    'linewidth',1.5,...
                                    'MarkerFaceColor',cols{i},...
                                    'MarkerFaceColor',cols{i},...
                                    'MarkerSize',7); hold on;
            end
            plot([min(y) max(y)],[0 0],'k--','linewidth',1.5);
            set(gca,'fontsize',25);
            set(gca,'linewidth',1.5);
            xlim([min(y) max(y)]);
            ylim([-0.1 0.02]);
            xlabel('$y_1/\sigma$','Interpreter','Latex','fontsize',25);
            ylabel('$\Pi \sigma^3/\varepsilon$','Interpreter','Latex','fontsize',25);
        end
        
        
        SaveFigure(['DisjoiningPressures'],res{k0}(1).config);        
    end    
    function PlotDensityProfiles(resC,config)
        CLT = ContactLineHS(config);
        CLT.PreprocessIDC();         

        figure('color','white','Position',[0 0 500 2000]); 
        subplot(3,1,1);
        CLT.IDC.do1DPlotNormal(resC{1}.rho1D,'s','b'); hold on;
        CLT.IDC.do1DPlotNormal(resC{2}.rho1D,'p','b');
        xlim([0 6]);
        xlabel('$y_2/\sigma$','Interpreter','Latex','fontsize',20);
        ylabel('$n \sigma^3$','Interpreter','Latex','fontsize',20);            

        subplot(3,1,2);
        for k = 3:length(resC)
            CLT.IDC.do1DPlotNormal(resC{k}.rho1D_WL,syms{k-2},'k'); hold on;
        end
        xlim([0 6]);
        ylim([0 1.8]);
        xlabel('$y_2/\sigma$','Interpreter','Latex','fontsize',20);
        ylabel('$n \sigma^3$','Interpreter','Latex','fontsize',20);    

        subplot(3,1,3);
        for k = 3:length(resC)
            CLT.IDC.do1DPlotNormal(resC{k}.rho1D_WG,syms{k-2},'k'); hold on;
        end
        xlim([0 6]);
        ylim([0 0.25]);
        xlabel('$y_2/\sigma$','Interpreter','Latex','fontsize',20);
        ylabel('$n \sigma^3$','Interpreter','Latex','fontsize',20);
        %SaveFigure('DensityProfiles',v2struct(config)); %N
    end   
    function res = PlotDensityProfile(res)
        
        CL = ContactLineHS(res.config);
        CL.Preprocess();
        CL.ComputeEquilibrium();                              
        CL.InitInterpolation([-12 12],[0.5 24.5]); close all;
        figure('Color','white','Position',[0 0 800 800]);
        CL.PlotContourResults({})
        set(gca,'fontsize',35);  
        xlabel('${y_1}/{\sigma}$','Interpreter','Latex','fontsize',35);
        ylabel('${y_2}/{\sigma}$','Interpreter','Latex','fontsize',35);
        CL.SaveCurrentFigure('EquilibriumContour');        
    end
    function res = ComputeError(in,h)
        conf = in.config;
        n    = in.N;                          
                
        for i = 1:length(n)
            
            conf.optsNum.PhysArea.N       = [n(i),n(i)];
            conf.optsNum.PhysArea.N2bound = max(10,2*round(n(i)/6));
                
            CL            = ContactLineHS(conf);
            preErr        = CL.Preprocess(); 
            
            res(i).config = conf;            

            %res(i).error_conv1 = preErr.error_conv1;
            %res(i).Conv        = CL.IntMatrV2.Conv;

           % [~,res(i).rho_1D_WL,params] = CL.Compute1D('WL');
%             res(i).error_wl = params.contactDensity_relError;
% 
%             [~,res(i).rho_1D_WG,params] = CL.Compute1D('WG');
%             res(i).error_wg = params.contactDensity_relError;
            
            CL.ComputeEquilibrium();      
            [res(i).sumRuleII_relError,...
             res(i).sumRuleII_eps,res(i).piII]  = CL.SumRuleIIError(in.Interval_SumRule);   
            %CLT.PostProcess(opts);
            %CLT.PlotDensitySlices();
            %CLT.PlotDisjoiningPressures();       

            close all;
            clear('CL');                            
        end
%        res.error_wl = error_wl;
%        res.error_wg = error_wg; 
    end     
    function PlotErrorGraph(resA,var_name,col,name)
        
        
                
        for k0 = 1:length(resA)
            n    = 1;
            resK = resA{k0};
            for k1 = 1:size(resK,1)
                for k2 = 1:size(resK,2)
                    line(n)   = abs(resK(k1,k2).(var_name));
                    line_N(n) = (resK(k1,k2).config.optsNum.PhysArea.N(1));%+res(k1,k2+1).NS)/2;
                    %plot(line_N(n),line(n),...
                     %       [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;
                    n = n+1;

                end
            end
            plot(line_N,line,...
                        ['-',syms{k0},col],'MarkerSize',10,'MarkerFaceColor',col); hold on;        
            %legendstring(end+1) = {name};
        end
                            
        
    end

end