function CheckSumRule_BH_HalfSpace()
    
    AddPaths('CodePaper');  
    close all;
    
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

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.);%1.375);%1.25)s;
    

    optsPhys = struct('V1',V1,...%'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);

    config = v2struct(optsNum,optsPhys);                        

    N    = 20:5:120;
    NS   = 80;%10:10:40;        
    
    ignoreList = {'config_optsNum_PhysArea_N',...
                  'config_optsNum_PhysArea_N2bound',...
                  'config_optsNum_V2Num_N',...
                  'config_optsNum_FexNum_N1disc',...
                  'config_optsNum_FexNum_N2disc'};
              
    comp = [];
    
    %**********************
    %**********************
    eta = 0.3;    
    res{1} = DataStorage('SumRuleError',@ComputeError,v2struct(N,NS,config,eta),[],comp,ignoreList);
    resC{1}.config = config;
    resC{1}.eta    = eta;
        	    
    eta            = 0.15;    
    res{2}         = DataStorage('SumRuleError',@ComputeError,v2struct(N,NS,config,eta),[],comp,ignoreList);  
    resC{2}.config = config;
    resC{2}.eta    = eta;
            
    config.optsNum.V2Num = struct('Fex','SplitAnnulus','N',[80,80]);
    config.optsPhys.V2   = struct('V2DV2','BarkerHendersonCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',5);
    config.optsPhys.V1.epsilon_w = 0.94;    
    res{3}        = DataStorage('SumRuleError',@ComputeError,v2struct(N,NS,config),[],comp,ignoreList); 
    resC{3}.config = config;
        
    config.optsNum.V2Num = struct('Fex','SplitAnnulus','N',[80,80]);
    config.optsPhys.V2   = struct('V2DV2','BarkerHendersonHardCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',5);   
    config.optsPhys.V1.epsilon_w = 1.0;    
    res{4} = DataStorage('SumRuleError',@ComputeError,v2struct(N,NS,config),[],comp,ignoreList); 
    resC{4}.config = config;        
    
    % **** 4 ****
    config.optsNum.V2Num  = struct('Fex','SplitDisk','N',[80,80]); 
    config.optsPhys.V2    = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 
    config.optsPhys.V1.epsilon_w = 1.0;    
    res{5} = DataStorage('SumRuleError',@ComputeError,v2struct(N,NS,config),[],comp,ignoreList);            
    resC{5}.config = config;       
    
    config.optsNum.V2Num     = struct('Fex','SplitDisk','N',[80,80]); 
    config.optsPhys.V2       = struct('V2DV2','ExponentialDouble','epsilon',1,'LJsigma',1);
    config.optsPhys.V1.V1DV1 = 'Vext_Exp_HardWall';
    config.optsPhys.V1.epsilon_w = 1.45;    
    res{6} = DataStorage('SumRuleError',@ComputeError,v2struct(N,NS,config),[],[],ignoreList);    
    resC{6}.config = config;    
    
    %**********************
    %**********************
    cols = GetGreyShades(length(N)); %cols = {'b','m','k','r','g'};
    syms = {'o','>','*','^','h'};        
    resC = DataStorage([],@ComputeDensityProfile,config,resC,[],ignoreList);
    
    PlotConvergenceDensityProfiles(res{1},N,[7 7 7]*1e-3,{[60,80],[80,100],[80 120]},'HS');    
    PlotConvergenceDensityProfiles(res{3},N,[5 5 5]*1e-2,{[60,80],[80,100],[100 120]},'BH');
    PlotConvergenceDensityProfiles(res{6},N,[7 7 7]*1e-2,{[60,80],[80,100],[100 120]},'exp');
    
    PlotDensityProfiles(resC);
         
    %**********************
    %**********************
    figure('color','white','Position',[0 0 900 800]); 
    legendstring = {};
    PlotErrorGraph(res{1},'error_wl','s-','b',['Hard sphere, eta = ',num2str(0.3)]); 
    PlotErrorGraph(res{2},'error_wl','p-','b',['Hard sphere, eta = ',num2str(0.15)]);     
    
    PlotErrorGraph(res{3},'error_wl','o-','k','BarkerHendersonCutoff_2D, liq');
    PlotErrorGraph(res{3},'error_wg','o--','k','BarkerHendersonCutoff_2D, vap');
    
    PlotErrorGraph(res{4},'error_wl','>-','k','BarkerHendersonHardCutoff_2D, liq');
    PlotErrorGraph(res{4},'error_wg','>--','k','BarkerHendersonHardCutoff_2D, vap');
    
    PlotErrorGraph(res{5},'error_wl','*-','k','BarkerHenderson_2D, liq');
    PlotErrorGraph(res{5},'error_wg','*--','k','BarkerHenderson_2D, vap');
        
    PlotErrorGraph(res{6},'error_wl','^-','k','ExponentialDouble, liq');
    PlotErrorGraph(res{6},'error_wg','^--','k','ExponentialDouble, vap');
    
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);
    xlabel('$N$','Interpreter','Latex','fontsize',20);
    ylabel(['Relative sum rule error $\frac{n(0)-n_C}{n_C}$'],...
                            'Interpreter','Latex','fontsize',20);
    xlim([(N(1)-2),(N(end)+2)]);    
    %legend(legendstring,'Location','northeast');%,'Orientation','horizontal');
   % legend(legendstring,'Location','eastOutside');%,'Orientation','horizontal');
    
    comment = 'Computed for hard wall, hard sphere fluid and BH fluid';    
    SaveFigure('SumRuleError',v2struct(N,NS,config,comment));
    
    function PlotConvergenceDensityProfiles(res,N,acc,ints,name)
        figure('Position',[0 0 1800 500],'color','white');
        subplot(1,3,1); PlotConvergenceDensityProfilesInt(res,N,ints{1},acc(1));
        subplot(1,3,2); PlotConvergenceDensityProfilesInt(res,N,ints{2},acc(2));
        subplot(1,3,3); PlotConvergenceDensityProfilesInt(res,N,ints{3},acc(3));
        saveC   = res(1,1).config;
        saveC.N = N;    
        SaveFigure(['DensityProfiles_Convergence_',name],saveC);
    end

    function PlotConvergenceDensityProfilesInt(res,N,NInt,yLimMax)
        conf = res.config;
        CLT = ContactLineHS(conf);
        CLT.PreprocessIDC();         

        xIP = (-1:0.002:1)';
        yIP = CLT.IDC.PhysSpace2(xIP);
        
        if(isfield(res(1,1),'rho_1D'))
            rhoName = 'rho_1D';
        else
            rhoName = 'rho_1D_WL';
        end
        
        rhoRef = res(1,1).(rhoName)(end);                
        mark   = ( (N >= NInt(1)) & (N <= NInt(2)));        
        is     = 1:size(res,1);
        is     = is(mark);        
        cls = GetGreyShades(length(is)); %cols = {'b','m','k','r','g'};
                
        for i = 1:sum(mark)
            resI   = res(is(i),1);
            rhoRefErr = abs(rhoRef-res(i,1).(rhoName)(end));
            if(rhoRefErr > 1e-13)
                error('reference density not unique');
            end
            disp(resI.config.optsNum.PhysArea.N(2))
            plot(resI.y2-0.5,resI.(rhoName)-rhoRef,'o','color',cls{i},'MarkerSize',8,'MarkerFaceColor',cls{i}); hold on;
            IP = barychebevalMatrix(resI.x2,xIP);
            plot(yIP-0.5,IP*resI.(rhoName)-rhoRef,'-','color',cls{i},'linewidth',1.5); hold on;
        end
        xlim([3 7]);
        ylim([-1 1]*yLimMax);
        set(gca,'fontsize',20);
        set(gca,'linewidth',1.5);
        xlabel('$y_2/\sigma$','Interpreter','Latex','fontsize',25);
        ylabel('$(n - n_{b})\sigma^3$','Interpreter','Latex','fontsize',25); 
       
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
        N      = in.N;
        NS     = in.NS;
        
        %error_wg    = zeros(length(N),length(NS));
        %error_wl    = zeros(length(N),length(NS));
                
        for i = 1:length(N)
            
            conf.optsNum.PhysArea.N       = [1,N(i)];
            conf.optsNum.PhysArea.N2bound = max(10,2*round(N(i)/6));
            
            for j = 1:length(NS)
                if(isfield(conf.optsPhys,'V2'))
                    conf.optsNum.V2Num.N       = [NS(j),NS(j)];
                end
                conf.optsNum.FexNum.N1disc = NS(j);
                conf.optsNum.FexNum.N2disc = NS(j);

                CL = ContactLineHS(conf);
                preErr = CL.Preprocess(); 
                
                res(i,j).config = conf;
                res(i,j).y2      = CL.IDC.Pts.y2;
                res(i,j).x2      = CL.IDC.Pts.x2;
                
                if(~isfield(conf.optsPhys,'V2'))                                        
                    [~,res(i,j).rho_1D,params] = CL.Compute1D(in.eta);
                    res(i,j).error_wl = params.contactDensity_relError;
                else
                    res(i,j).error_conv1 = preErr.error_conv1;
                    res(i,j).Conv        = CL.IntMatrV2.Conv;
                    
                    [~,res(i,j).rho_1D_WL,params] = CL.Compute1D('WL');
                    res(i,j).error_wl = params.contactDensity_relError;

                    [~,res(i,j).rho_1D_WG,params] = CL.Compute1D('WG');
                    res(i,j).error_wg = params.contactDensity_relError;
                end
                %res(i,j).error_n2_1    = preErr.error_n2_1;
                %res(i,j).error_n3_1    = preErr.error_n3_1;
                %res(i,j).error_n2v2_1  = preErr.error_n2v2_1;                                
                
                res(i,j).N  = N(i);
                res(i,j).NS = NS(j); 
                %err         = CheckMeanfieldConvolution(CL);
                %%res(i,j).error_conv1 = err.error_conv1;                 
                
                %res(i,j).A_n2       = CL.IntMatrFex.AD.n2;
                %res(i,j).A_n3       = CL.IntMatrFex.AD.n3;
                %res(i,j).A_n2_v_1   = CL.IntMatrFex.AD.n2_v_1;
                %res(i,j).A_n2_v_2   = CL.IntMatrFex.AD.n2_v_2;                                
                
                close all;
                clear('CL');                
            end
        end
%        res.error_wl = error_wl;
%        res.error_wg = error_wg; 
    end     
    function PlotErrorGraph(res,var_name,sym,col,name)
        
        n = 1;
                
        for k1 = 1:size(res,1)
            for k2 = 1:size(res,2)
                line(n)   = abs(res(k1,k2).(var_name));
                line_N(n) = (res(k1,k2).N);%+res(k1,k2+1).NS)/2;
                %plot(line_N(n),line(n),...
                 %       [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;
                n = n+1;
                
                                
            end
        end
        plot(line_N,line,...
                        [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;        
                    
        legendstring(end+1) = {name};
    end

end