function ThesisNumerics_CheckSumRule_BH_HalfSpace()
    
    AddPaths('CodePaper');  
    close all;
    global loadAll
    loadAll = true;
    %cols = GetGreyShades(length(N)); %cols = {'b','m','k','r','g'};
    syms = {'s','p','o','^','>','*','h'};
    %lines = {'-','--'};
    
    PhysArea = struct('N',[1,40],...
                      'L1',5,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',90);
                      
    %V2Num   = struct('Fex','SplitDisk','N',[34,34]);
    %V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 
    %    V2 = struct('V2DV2','Exponential','epsilon',1.5,'LJsigma',1);     

    FexNum   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',80,'N2disc',80);                   

    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',FexNum,...%'V2Num',V2Num,...
                     'maxComp_y2',-1,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.);%1.375);%1.25);

    optsPhys = struct('V1',V1,...%'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);

    config = v2struct(optsNum,optsPhys);                        

    N    = 20:5:120;        
    y    = (0.5:0.1:20)';
    
    ignoreList = {'config_optsNum_PhysArea_N',...
                  'config_optsNum_PhysArea_N2bound',...
                  'config_optsNum_V2Num_N',...
                  'config_optsNum_FexNum_N1disc',...
                  'config_optsNum_FexNum_N2disc',...
                  'NS'};
              
    comp = [];
    
    %**********************
    %**********************
    %eta = 0.3;
    eta = 0.3257;
    res{1} = DataStorage('SumRuleError',@ComputeError,v2struct(N,config,eta),[],comp,ignoreList);    
    resC{1}.config = config;
    resC{1}.name = 'HSeta0.3';
    resC{1}.eta    = eta;
        	    
    %eta            = 0.15;    
    eta            = 0.0147;    
    res{2}         = DataStorage('SumRuleError',@ComputeError,v2struct(N,config,eta),[],comp,ignoreList);      
    resC{2}.config = config;
    resC{2}.name = 'HSeta0.15';
    resC{2}.eta    = eta;    
                
    config.optsNum.V2Num = struct('Fex','SplitAnnulus','N',[80,80]);
    %config.optsNum.V2Num = struct('Fex','SplitDisk','N',[80,80]);
    %config.optsPhys.V2   = struct('V2DV2','BarkerHendersonCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',5.0);
    config.optsPhys.V2   = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);%2.5
    config.optsPhys.V1.epsilon_w = 0.865; %0.928;%94;% 0.94;    	
    res{3}        = DataStorage('SumRuleError',@ComputeError,v2struct(N,config),[],[],ignoreList); 
    resC{3}.name = 'BH';
    resC{3}.config = config;           
    
%     config.optsNum.V2Num     = struct('Fex','SplitDisk','N',[80,80]); 
%     config.optsPhys.V2       = struct('V2DV2','ExponentialDouble','epsilon',1,'LJsigma',1);
%     config.optsPhys.V1.V1DV1 = 'Vext_Exp_HardWall';
%     config.optsPhys.V1.epsilon_w = 1.45;    
%     res{4} = DataStorage('SumRuleError',@ComputeError,v2struct(N,config),[],[],ignoreList);    
%     resC{4}.name = 'exp';
%     resC{4}.config = config;               
    
    %**********************
    %**********************
    res = PostProcess(res,y);    
    
 %   resC = DataStorage([],@ComputeDensityProfile,config,resC,true,ignoreList);
             
    %**********************
    %**********************        
    %PlotConvergenceDensityProfilesInt(res{3},N,0.17,[1 7],[60,100],'rho_1D_WL'); title([]);
    PlotConvergenceDensityProfilesInt(res{3},N,0.05*[-1 1],[3 7],[60,100],'rho_1D_WL'); title([]);
    SaveFigure('BH_ConvergingDensityProfiles');
        
    PlotConvergenceDensityProfilesInt(res{1},N,[-0.0075 0.015],[3 7],[60,100],'rho_1D'); title([]);    
    SaveFigure('HS_ConvergingDensityProfiles');
    
%     PlotErrorVsY_Full(res,resC);
%     
%     PlotFullErrorGraph(res,'adsorption_err_wl','Relative adsorption error','AdsorptionError');
%     PlotFullErrorGraph(res,'error_wl',...
%                 'Relative sum rule error $\frac{n(0)-n_C}{n_C}$',...
%                 'SumRuleError');        
%     
%         
%     PlotFullErrorGraph(res,'rho_y_err_max_wl','$\max_y |n_N(y) - n_{N+\Delta N}(y)|$','MaxDensityError_Y');
%     %**********************
%     %**********************
%     PlotConvergenceDensityProfiles(res{3},N,[],[0 100],{[60,80],[80,100],[100 120]},'BH_WG','rho_1D_WG');
%     
%     PlotConvergenceDensityProfiles(res{1},N,[7 7 7]*1e-3,[3 7],{[60,80],[80,100],[80 120]},'HS','rho_1D');    
%     PlotConvergenceDensityProfiles(res{3},N,[5 5 5]*1e-2,[3 7],{[60,80],[80,100],[100 120]},'BH_WL','rho_1D_WL');

    
    function res = PostProcess(res,yIP)
        
        CLT = ContactLineHS(res{1}(1,1).config);
        CLT.PreprocessIDC(); 
        
        for i0 = 1:length(res)
            for i1 = 1:size(res{i0},1)
                for i2 = 1:size(res{i0},2)
                    resI = res{i0}(i1,i2);
                    
                    [x,w]  = ClenCurtFlip(length(resI.x2)-1);
                    [yy,dy] = CLT.IDC.PhysSpace2(x);
                    Int2   = dy'.*w;
                    Int2(Int2==inf) = 0;
                    %Int2(yy>) = 0;

                    xx   = CLT.IDC.CompSpace2(yIP);
                    IP   = barychebevalMatrix(resI.x2,xx);
                    res{i0}(i1,i2).y_IP = yIP;
                    res{i0}(i1,i2).Int2 = Int2;
                    if(isfield(resI,'rho_1D'))
                        res{i0}(i1,i2).rho_y_wl      = IP*resI.rho_1D;
                        res{i0}(i1,i2).adsorption_wl = Int2*(resI.rho_1D-resI.rho_1D(end));
                    end
                    if(isfield(resI,'rho_1D_WL'))
                        res{i0}(i1,i2).rho_y_wl      = IP*resI.rho_1D_WL;
                        res{i0}(i1,i2).adsorption_wl = Int2*(resI.rho_1D_WL-resI.rho_1D_WL(end));
                    end
                    if(isfield(resI,'rho_1D_WG'))
                        res{i0}(i1,i2).rho_y_wg      = IP*resI.rho_1D_WG;
                        res{i0}(i1,i2).adsorption_wg = Int2*(resI.rho_1D_WG-resI.rho_1D_WG(end));
                    end
                end
            end
        end
                
        for i0 = 1:length(res)
            for i1 = 1:(length(res{i0})-1)  
                res{i0}(i1).adsorption_err_wl = abs(res{i0}(i1).adsorption_wl-res{i0}(i1+1).adsorption_wl)/res{i0}(i1+1).adsorption_wl;
                res{i0}(i1).rho_y_err_wl      = abs(res{i0}(i1).rho_y_wl-res{i0}(i1+1).rho_y_wl);
                res{i0}(i1).rho_y_err_max_wl  = max(abs(res{i0}(i1).rho_y_err_wl));               

                if(isfield(res{i0}(i1),'rho_y_wg'))
                    res{i0}(i1).adsorption_err_wg = abs(res{i0}(i1).adsorption_wg - res{i0}(i1+1).adsorption_wg)/res{i0}(i1+1).adsorption_wg;
                    res{i0}(i1).rho_y_err_wg      = abs(res{i0}(i1).rho_y_wg-res{i0}(i1+1).rho_y_wg);
                    res{i0}(i1).rho_y_err_max_wg  = max(abs(res{i0}(i1).rho_y_err_wg));                        
                end                
            end
        end
    end    
    function PlotConvergenceDensityProfiles(res,N,acc,xlims,ints,name,varName)        
        figure('Position',[0 0 1800 500],'color','white');
        for i = 1:3            
            subplot(1,3,i); 
            if(~isempty(acc))
                PlotConvergenceDensityProfilesInt(res,N,acc(i),xlims,ints{i},varName);
            else
                PlotConvergenceDensityProfilesInt(res,N,[],xlims,ints{i},varName);
            end
        end        
        saveC   = res(1,1).config;
        saveC.N = N;    
        SaveFigure(['DensityProfiles_Convergence_',name],saveC);
    end
    function PlotConvergenceDensityProfilesInt(res,N,yLimMax,xlims,NInt,rhoName)
        yMax = 10;
        
        conf = res.config;
        CLT = ContactLineHS(conf);
        CLT.PreprocessIDC();         

        yIP = 0.5+ (0:0.02:yMax)';
        xIP = CLT.IDC.CompSpace2(yIP);
        %xIP = (-1:0.002:1)';
        %yIP = CLT.IDC.PhysSpace2(xIP);  
        
        %mark_Y                = (yIP <= yMax + 0.5);
        %mark_Y(sum(mark_Y)+1) = true;
        %xIP                   = xIP(mark_Y);
        %yIP                   = yIP(mark_Y);
        
        rhoRef = res(1).(rhoName)(end);                
        mark   = ( (N >= NInt(1)) & (N <= NInt(2)));        
        is     = 1:length(res);
        is     = is(mark);        
        cls = GetGreyShades(length(is)); %cols = {'b','m','k','r','g'};
                
        f1 = figure('color','white','Position',[0 0 250 200]);
       %plot(res(end).y2-0.5,res(end).(rhoName),'o','color','k','MarkerSize',8,'MarkerFaceColor','k'); hold on;
        IP = barychebevalMatrix(res(end).x2,xIP);
        [xs,ys] = GetReasonableRange(yIP-0.5,IP*res(end).(rhoName));
        plot(xs,ys,'-','color','k'); hold on;
        xlim([0 10])        
        xlabel('$y_2$','Interpreter','Latex'); %/\sigma
        ylabel('$n$','Interpreter','Latex'); %\sigma^3
        
        f2 = figure('color','white');
        for i = 1:sum(mark)
            resI   = res(is(i));
            rhoRefErr = abs(rhoRef-res(i).(rhoName)(end));
            if(rhoRefErr > 1e-13)
                error('reference density not unique');
            end
            %disp(resI.config.optsNum.PhysArea.N(2))
            mark_Y2 = (resI.y2 <= yMax + 0.5);
            rho_o   = resI.(rhoName);
            [xs,ys] = GetReasonableRange(resI.y2(mark_Y2)-0.5,rho_o(mark_Y2)-rhoRef,yLimMax);
            plot(xs,ys,'o','color',cls{i},'MarkerFaceColor',cls{i}); hold on;
            IP = barychebevalMatrix(resI.x2,xIP);
            [xs,ys] = GetReasonableRange(yIP-0.5,IP*resI.(rhoName)-rhoRef,yLimMax);
            plot(xs,ys,'-','color',cls{i}); hold on;
        end        
        xlim(xlims);
        if(~isempty(yLimMax))                    
            ylim(yLimMax);
        else
            set(gca,'YScale','log');            
            ylim([1e-12 2]);
        end
      %  title(['$N \in [',num2str(NInt(1)),';',num2str(NInt(2)),']$'],...
%                                    'Interpreter','Latex',...
%                                    'fontsize',20);        
        xlabel('$y_2$','Interpreter','Latex'); %/\sigma
        ylabel('$n - n_{b}$','Interpreter','Latex'); %\sigma^3        
       
        inset2(f1,f2,0.5,[0.5,0.5]);   close(f2);    
    end    
    function PlotDensityProfiles(resC)
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
        ylim([0 3]);
        xlabel('$y_2/\sigma$','Interpreter','Latex','fontsize',20);
        ylabel('$n \sigma^3$','Interpreter','Latex','fontsize',20);    

        subplot(3,1,3);
        for k = 3:length(resC)
            CLT.IDC.do1DPlotNormal(resC{k}.rho1D_WG,syms{k-2},'k'); hold on;
        end
        xlim([0 6]);
        ylim([0 0.6]);
        xlabel('$y_2/\sigma$','Interpreter','Latex','fontsize',20);
        ylabel('$n \sigma^3$','Interpreter','Latex','fontsize',20);
        SaveFigure('DensityProfiles',v2struct(N,NS,config));
    end   
    function res = ComputeDensityProfile(conf,res)
        
        CL = ContactLineHS(res{1}.config);
        CL.Preprocess(); 
        
        for i = 1:length(res)
            
            if(isfield(res{i}.config.optsNum,'V2Num'))                
                CL.optsPhys.V1.V1DV1     =  res{i}.config.optsPhys.V1.V1DV1;
                CL.optsPhys.V1.epsilon_w = res{i}.config.optsPhys.V1.epsilon_w; 
                CL.optsNum.V2Num         = res{i}.config.optsNum.V2Num;
                CL.optsPhys.V2           = res{i}.config.optsPhys.V2;            
                CL.Preprocess_MeanfieldContribution();
                
                [CL.optsPhys.rhoGas_sat,...
                 CL.optsPhys.rhoLiq_sat,...
                 CL.optsPhys.mu_sat,this.optsPhys.p] = BulkSatValues(CL.optsPhys);
                
                [~,res{i}.rho1D_WL] = CL.Compute1D('WL');
                [~,res{i}.rho1D_WG] = CL.Compute1D('WG');
            else
                [~,res{i}.rho1D] = CL.Compute1D(res{i}.eta);
            end          
            close all;
        end
        close all;
    end
    function res = ComputeError(in,h)
        conf = in.config;
        n    = in.N;                        
                
        for i = 1:length(n)
            
            conf.optsNum.PhysArea.N       = [1,n(i)];
            conf.optsNum.PhysArea.N2bound = max(10,2*round(n(i)/6));

            CL = ContactLineHS(conf);
            preErr = CL.Preprocess(); 

            [~,~,res(i).Int_1D]  = CL.IDC.ComputeIntegrationVector();

            res(i).config = conf;
            res(i).y2      = CL.IDC.Pts.y2;
            res(i).x2      = CL.IDC.Pts.x2;

            if(~isfield(conf.optsPhys,'V2'))                                        
                [~,res(i).rho_1D,params] = CL.Compute1D(in.eta);
                res(i).error_wl = params.contactDensity_relError;                                                                                
            else
                res(i).error_conv1 = preErr.error_conv1;
                res(i).Conv        = CL.IntMatrV2.Conv;

                [~,res(i).rho_1D_WL,params] = CL.Compute1D('WL');
                res(i).error_wl = params.contactDensity_relError;

                [~,res(i).rho_1D_WG,params] = CL.Compute1D('WG');
                res(i).error_wg = params.contactDensity_relError;
            end

            res(i).N  = n(i);                
                
            close all;
            clear('CL');                            
        end
%        res.error_wl = error_wl;
%        res.error_wg = error_wg; 
    end     

    function PlotErrorVsY_Full(res,resC)
        figure('color','white','Position',[0 0 1400 1400]); 
        subplot(2,2,1); PlotErrorVsY(res{1},'rho_y_err_wl','-o',resC{1}.name);
        subplot(2,2,2); PlotErrorVsY(res{2},'rho_y_err_wl','-o',resC{2}.name);
        subplot(2,2,3); PlotErrorVsY(res{3},'rho_y_err_wl','-o',[resC{3}.name,'_{wl}']);
        subplot(2,2,4); PlotErrorVsY(res{3},'rho_y_err_wg','-o',[resC{3}.name,'_{wg}']);
%        subplot(3,2,5); PlotErrorVsY(res{4},'rho_y_err_wl','-o',[resC{4}.name,'_{wl}']);
%        subplot(3,2,6); PlotErrorVsY(res{4},'rho_y_err_wg','-o',[resC{4}.name,'_{wg}']);
        SaveFigure('ConvergenceError_Y');
    end
    function PlotErrorVsY(res,var_name,sym,title_name)        
        ns   = 1:length(res);
        mark = 5:9;%length(res);
        ns   = ns(mark);
        cls  = GetGreyShades(length(ns));
        for j1 = 1:length(ns)                                    
            col = cls{j1};
            k1  = ns(j1); 
            disp(res(k1).N);
            if(isfield(res(k1),var_name) && ~isempty(res(k1).(var_name)))
                plot(res(k1).y_IP,res(k1).(var_name),sym,'color',col,...
                        'MarkerSize',10,'MarkerFaceColor',col);
                hold on;
            else
                break;
            end            
        end        
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$y_2$','Interpreter','Latex','fontsize',20);
        title(title_name,'fontsize',20);
    end
       
    function PlotFullErrorGraph(res,var_name,yLab,filename)
        figure('color','white','Position',[0 0 900 800]); 
        for i0 = 1:length(res)            
            var_name_wg = [var_name(1:end-3),'_wg'];
            if(isfield(res{i0}(1,1),var_name_wg))
                PlotErrorGraph(res{i0},var_name,['-',syms{i0}],'k');
                PlotErrorGraph(res{i0},var_name_wg,['--',syms{i0}],'k');
            else
                PlotErrorGraph(res{i0},var_name,['-',syms{i0}],'k');                
            end
        end
        set(gca,'YScale','log');
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);
        xlabel('$N$','Interpreter','Latex','fontsize',20);
        ylabel(yLab,'Interpreter','Latex','fontsize',20);
        xlim([(N(1)-2),(N(end)+2)]);    
        SaveFigure(filename);
    end
    function PlotErrorGraph(res,var_name,sym,col)
        
        n = 1;                
        for k1 = 1:length(res)
            if(isfield(res(k1),var_name) && ~isempty(res(k1).(var_name)))
                line(n)   = abs(res(k1).(var_name));
                line_n(n) = (res(k1).N);%+res(k1,k2+1).NS)/2;                
                n = n+1;
            else
                break;
            end
        end
        plot(line_n,line,...
                        [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;        
    end
    function [x,y] = GetReasonableRange(x,y,yLimMax)
        
        yM             = max(abs(y));
        y(abs(y) < yM*1e-4) = 0;       
        
        if(nargin >= 3)
            mark = ((y >= yLimMax(1)) & (y <= yLimMax(2)));
        else
            mark = true(size(x));
        end
        x = x(mark);
        y = y(mark);
    end

end