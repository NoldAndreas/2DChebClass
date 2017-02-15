function JCP_639_2017_Fig7_8_CheckSumRule()
    
    AddPaths('JCP_639_2017');  
    close all;
    
    %*************************
    % **** Initialization ****
    %*************************            
    PhysArea = struct('N',[1,40],...
                      'L1',5,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',90);                      
    FexNum   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',80,'N2disc',80);                   
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',FexNum,...
                     'maxComp_y2',-1,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.);
    optsPhys = struct('V1',V1,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);

    config = v2struct(optsNum,optsPhys);                        
           
    %**********************************
    % **** Computations for Fig. 7 ****
    %**********************************
    nn  = 50:25:150;
    res = DoComputations(nn,config);
    
    PlotConvergenceDensityProfilesInt(res{3},nn,0.05*[-1 1],[3 7],[50,150],'rho_1D_WL'); title([]);    
    SaveFigure('JCP_639_2017_Fig7_a');
        
    PlotConvergenceDensityProfilesInt(res{1},nn,[-0.0075 0.015],[3 7],[50,150],'rho_1D'); title([]);    
    SaveFigure('JCP_639_2017_Fig7_b');
    
    %**********************************
    % **** Computations for Fig. 8 ****
    %**********************************    
    res = DoComputations(50:50:250,config);
           
    PlotFullErrorGraph(res,'error_wl','Relative sum rule error $\frac{n(0)-n_C}{n_C}$');        
    SaveFigure('JCP_639_2017_Fig8_a');
            
    PlotFullErrorGraph(res,'rho_y_err_max_wl','$\max_y |n_N(y) - n_{N+\Delta N}(y)|$');
    SaveFigure('JCP_639_2017_Fig8_b');
    
    function res = DoComputations(N,config)
           ignoreList = {'config_optsNum_PhysArea_N',...
                  'config_optsNum_PhysArea_N2bound',...
                  'config_optsNum_V2Num_N',...
                  'config_optsNum_FexNum_N1disc',...
                  'config_optsNum_FexNum_N2disc',...
                  'NS'};              
        
        eta            = 0.3257;
        res{1}         = DataStorage('SumRuleError',@ComputeError,v2struct(N,config,eta),[],[],ignoreList);    

        eta            = 0.0147;    
        res{2}         = DataStorage('SumRuleError',@ComputeError,v2struct(N,config,eta),[],[],ignoreList);      

        config.optsNum.V2Num = struct('Fex','SplitAnnulus','N',[80,80]);
        config.optsPhys.V2   = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);
        config.optsPhys.V1.epsilon_w = 0.865;
        res{3}        = DataStorage('SumRuleError',@ComputeError,v2struct(N,config),[],[],ignoreList);                       

        y   = (0.5:0.1:20)';
        res = PostProcess(res,y);    
    end
    
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
    function PlotConvergenceDensityProfilesInt(res,N,yLimMax,xlims,NInt,rhoName)
        conf = res.config;
        CLT = ContactLineHS(conf);
        CLT.PreprocessIDC();         

        yIP = (0.5:0.01:11)';
        xIP = CLT.IDC.CompSpace2(yIP);
        
        rhoRef = res(1).(rhoName)(end);                
        mark   = ( (N >= NInt(1)) & (N <= NInt(2)));        
        is     = 1:length(res);
        is     = is(mark);        
        cls = GetGreyShades(length(is)); 
                
        f1 = figure('color','white');
        IP = barychebevalMatrix(res(end).x2,xIP);
        plot(yIP-0.5,IP*res(end).(rhoName),'-','color','k','linewidth',1.5); hold on;
        xlim([0 10])
        set(gca,'fontsize',20); set(gca,'linewidth',1.5);
        xlabel('$y_2$','Interpreter','Latex','fontsize',25);
        ylabel('$n$','Interpreter','Latex','fontsize',25); 
        
        f2 = figure('color','white');
        for i = 1:sum(mark)
            resI   = res(is(i));
            rhoRefErr = abs(rhoRef-res(i).(rhoName)(end));
            if(rhoRefErr > 1e-13)
                error('reference density not unique');
            end
            IP = barychebevalMatrix(resI.x2,xIP);
            plot(yIP-0.5,IP*resI.(rhoName)-rhoRef,'-','color',cls{i},'linewidth',1.5); hold on;
        end        
        xlim(xlims);
        if(~isempty(yLimMax))                    
            ylim(yLimMax);
        else
            set(gca,'YScale','log');            
            ylim([1e-12 2]);
        end
        set(gca,'fontsize',15); set(gca,'linewidth',1.5);
        xlabel('$y_2$','Interpreter','Latex','fontsize',20);
        ylabel('$n - n_{b}$','Interpreter','Latex','fontsize',20);
       
        inset2(f1,f2,0.5,[0.5,0.5]);   close(f2);    
    end    
    function res = ComputeError(in,h)
        conf = in.config;
        n    = in.N;                        
                
        for i = 1:length(n)
            
            conf.optsNum.PhysArea.N       = [1,n(i)];
            conf.optsNum.PhysArea.N2bound = max(10,2*round(n(i)/6));
            conf.optsNum.FexNum.N1disc    =  2*round((20+n(i)/4)/2)*[1,1];
            conf.optsNum.FexNum.N2disc    =  2*round((20+n(i)/4)/2)*[1,1];
            if(isfield(conf.optsNum,'V2Num'))
                conf.optsNum.V2Num.N      =  2*round((20+n(i)/2)/2)*[1,1];
            end

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
    end            
    function PlotFullErrorGraph(res,var_name,yLab)
        syms = {'s','p','o','^','>','*','h'};        
        figure('color','white');
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
    end
    function PlotErrorGraph(res,var_name,sym,col)        
        n = 1;                
        for k1 = 1:length(res)
            if(isfield(res(k1),var_name) && ~isempty(res(k1).(var_name)))
                line(n)   = abs(res(k1).(var_name));
                line_n(n) = (res(k1).N);
                n = n+1;
            else
                break;
            end
        end
        plot(line_n,line,...
                        [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;        
    end

end