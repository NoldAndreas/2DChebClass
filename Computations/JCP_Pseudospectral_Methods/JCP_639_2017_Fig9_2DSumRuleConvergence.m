function JCP_639_2017_Fig9_2DSumRuleConvergence()
    
    AddPaths('JCP_639_2017');
    close all;
    
    %**************************
    % **** Initialization *****
    %**************************
    PhysArea = struct('N',[20,20],...
                      'L1',4,'L2',2,...
                      'alpha_deg',90);                      
    V2Num    = struct('Fex','SplitAnnulus','N',[80,80]);    
    FexNum   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);
    optsNum  = struct('PhysArea',PhysArea,...
                      'FexNum',FexNum,'V2Num',V2Num,...
                      'maxComp_y2',25,...
                      'y1Shift',0);

    V1       = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.856);
    V2       = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);         
    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);

    config   = v2struct(optsNum,optsPhys);                        

    Interval_SumRule = [-10,20];
    N                = 20:10:60;        
    
    ignoreList = {'config_optsNum_PhysArea_N',...
                  'config_optsNum_PhysArea_N2bound'};        
              
	%************************
    % **** Computations *****
    %************************            	                
    res{1} = DataStorage('SumRuleError2D',@ComputeError,v2struct(N,config,Interval_SumRule),[],[],ignoreList);    
    res    = PostProcess(res);
    
    %**************************
    % **** Plotting *****
    %**************************
    
    % **** Plot contour plot ****
    PlotDensityProfile(res{1}(3));
    
    % **** Plot disjoining pressure profiles ****
    PlotDisjoiningPressureProfiles(res);
    SaveFigure('JCP_639_2017_Fig9_a',v2struct(N,config));
    close all;
        
    % **** Plot computation times ****
    f0 = figure('color','white','Position',[0 0 800 800]); 
    PlotErrorGraph(res,'compTimeEq','k','k');
    set(gca,'yscale','log'); 
    ylim([50 10^4]);
    set(gca,'fontsize',15); set(gca,'linewidth',1.5);
    xlabel('$N$','Interpreter','Latex','fontsize',20);%/{\sigma}
	ylabel('t (sec.)','Interpreter','Latex','fontsize',20); %/{\sigma}    
    SaveFigure('CompTime');      
    
    % **** Plot error convergence ****
    f2 = figure('color','white','Position',[0 0 800 800]); 
    PlotErrorGraph(res,'sumRuleII_epsMax','-k','k');     
    PlotErrorGraph(res,'sumRuleII_relError','k');                
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);            
    xlabel('$N$','Interpreter','Latex','fontsize',20);
    xlim([(N(1)-2),(N(end)+2)]);    
	
    % **** Produce full figure ****
    inset2(f2,f0,0.3,[0.25,0.2]); close(f0);
    SaveFigure('JCP_639_2017_Fig9_b',v2struct(N,config));
    
    function PlotDisjoiningPressureProfiles(res)
        K    = length(res);
        y    = res{1}(1).piII.y1.Pts.y;
        cols = GetGreyShades(length(res{1}));
        syms = {'o','*','^','<','d','s','>'};      
        figure('color','white','Position',[0 0 800 800]);         
        for k0 = 1:K
            subplot(K,1,k0);
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
            set(gca,'fontsize',20);
            set(gca,'linewidth',1.5);
            xlim([min(y) max(y)]);
            ylim([-0.1 0.02]);
            xlabel('$y_1$','Interpreter','Latex','fontsize',20);
            ylabel('$\Pi$','Interpreter','Latex','fontsize',20);
        end        
    end    

    function PlotDensityProfiles(resC,config)
        CLT = ContactLineHS(config);
        CLT.PreprocessIDC();         

        figure('color','white','Position',[0 0 600 600]); 
        %subplot(3,1,1);
        CLT.IDC.do1DPlotNormal(resC{1}.rho1D,'s','b'); hold on;
        CLT.IDC.do1DPlotNormal(resC{2}.rho1D,'p','b');
        xlim([0 6]);
        xlabel('$y_2$','Interpreter','Latex','fontsize',20);
        ylabel('$n$','Interpreter','Latex','fontsize',20);  

        figure('color','white','Position',[0 0 600 600]);  
        for k = 3:length(resC)
            CLT.IDC.do1DPlotNormal(resC{k}.rho1D_WL,syms{k-2},'k'); hold on;
        end
        xlim([0 6]);
        ylim([0 1.8]);
        xlabel('$y_2$','Interpreter','Latex','fontsize',20); 
        ylabel('$n$','Interpreter','Latex','fontsize',20);   

        figure('color','white','Position',[0 0 600 600]);  %subplot(3,1,3);
        for k = 3:length(resC)
            CLT.IDC.do1DPlotNormal(resC{k}.rho1D_WG,syms{k-2},'k'); hold on;
        end
        xlim([0 6]);
        ylim([0 0.25]);
        xlabel('$y_2/\sigma$','Interpreter','Latex','fontsize',20);
        ylabel('$n \sigma^3$','Interpreter','Latex','fontsize',20);
        %SaveFigure('DensityProfiles',v2struct(config)); %N
    end   
    function fullname = PlotDensityProfile(res)
        CL = ContactLineHS(res.config);
        CL.Preprocess();
        CL.ComputeEquilibrium();
        LP = 15;
        CL.InitInterpolation([-1 1]*0.5*LP,[0 LP]+0.5); close all;
        figure('Color','white','Position',[0 0 800 800]);
        CL.PlotContourResults({})
        set(gca,'fontsize',35);  
        xlabel('${y_1}$','Interpreter','Latex','fontsize',35);
        ylabel('${y_2}$','Interpreter','Latex','fontsize',35);
        fullname = CL.SaveCurrentFigure('EquilibriumContour');                
    end
    function res = ComputeError(in,h)
        conf = in.config;
                
        for i = 1:length(in.N)
            
            conf.optsNum.PhysArea.N  = in.N(i)*[1,1];                
            CL            = ContactLineHS(conf);
            CL.Preprocess(); 
            
            res(i).config = conf;            
            r             = CL.ComputeEquilibrium(struct('solver','Picard'));      
            res(i).compTimeEq = r.compTime;
            [res(i).sumRuleII_relError,...
             res(i).sumRuleII_eps,res(i).piII]  = CL.SumRuleIIError(in.Interval_SumRule);   
         
            close all;
            clear('CL');                            
        end
    end     
    function PlotErrorGraph(resA,var_name,col,colMark)
        syms = {'o','*','^','<','d','s','>'};      
        if(nargin < 4)
            colMark = col;
        end
        for k0 = 1:length(resA)
            n    = 1;
            resK = resA{k0};
            for k1 = 1:size(resK,1)
                for k2 = 1:size(resK,2)
                    line(n)   = abs(resK(k1,k2).(var_name));
                    line_N(n) = (resK(k1,k2).config.optsNum.PhysArea.N(1));%+res(k1,k2+1).NS)/2;
                    n = n+1;
                end
            end
            plot(line_N,line,...
                        [syms{k0},'-',col],'MarkerSize',10,'MarkerFaceColor',colMark); 
            hold on;
        end        
    end

    function res = PostProcess(res)
        for i0 = 1:length(res)
            for i1 = 1:length(res{i0})
                DP = res{i0}(i1).piII.DP;
                res{i0}(i1).sumRuleII_epsMax = max(abs(DP(1)),abs(DP(end)))/max(abs(DP));
                
                if(ischar(res{i0}(i1).compTimeEq))
                    res{i0}(i1).compTimeEq = str2double(res{i0}(i1).compTimeEq);
                end
            end
        end
    end

end