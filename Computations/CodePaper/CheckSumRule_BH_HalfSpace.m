function CheckSumRule_BH_HalfSpace()
    
    AddPaths();    

    PhysArea = struct('N',[1,60],...
                      'L1',5,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',16,'h',1,...
                      'alpha_deg',90);
                      

    %V2Num   = struct('Fex','SplitDisk','N',[34,34]);
    %V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 
    %    V2 = struct('V2DV2','Exponential','epsilon',1.5,'LJsigma',1);     

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...%'V2Num',V2Num,...
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
    
    figure('color','white','Position',[0 0 900 800]); 
    legendstring = {};
    
    % **** 1 ****    
    eta = 0.3;    
    res{1} = DataStorage([],@ComputeError,v2struct(N,NS,config,eta),[]);
        
	AddPaths();    
    ChangeDirData();
    % **** 2 ****    
    eta = 0.15;    
    res{2} = DataStorage([],@ComputeError,v2struct(N,NS,config,eta),[]);  
    
    AddPaths();  ChangeDirData();   
    % **** 3 ****
   % NS   = 80;
   % N    = 20:10:120;
    config.optsNum.V2Num = struct('Fex','SplitAnnulus','N',[80,80]);
    config.optsPhys.V2   = struct('V2DV2','BarkerHendersonCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',5);
    config.optsPhys.V1.epsilon_w = 0.9;    
    res{3} = DataStorage([],@ComputeError,v2struct(N,NS,config),[]); 
    
    AddPaths();  ChangeDirData();     
    config.optsNum.V2Num = struct('Fex','SplitAnnulus','N',[80,80]);
    config.optsPhys.V2   = struct('V2DV2','BarkerHendersonHardCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',5);   
    config.optsPhys.V1.epsilon_w = 0.9;    
    res{4} = DataStorage([],@ComputeError,v2struct(N,NS,config),[]); 
    
    AddPaths();  ChangeDirData();
    
    % **** 4 ****
    config.optsNum.V2Num  = struct('Fex','SplitDisk','N',[80,80]); 
    config.optsPhys.V2    = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 
    config.optsPhys.V1.epsilon_w = 0.9;    
    res{5} = DataStorage([],@ComputeError,v2struct(N,NS,config),[]);            
    
    AddPaths();  ChangeDirData();
    
    config.optsNum.V2Num  = struct('Fex','SplitDisk','N',[80,80]); 
    config.optsPhys.V2    = struct('V2DV2','ExponentialDouble','epsilon',1,'LJsigma',1);
    config.optsPhys.V1.epsilon_w = 0.9;    
    res{6} = DataStorage([],@ComputeError,v2struct(N,NS,config),[]);    
    
    AddPaths();  ChangeDirData();
    
    
    PlotErrorGraph(res{1},'error_wl','o-','k',['Hard sphere, eta = ',num2str(0.3)]); 
    PlotErrorGraph(res{2},'error_wl','s-','k',['Hard sphere, eta = ',num2str(0.15)]);     
    
    PlotErrorGraph(res{3},'error_wl','d-','b','BarkerHendersonCutoff_2D, liq');
    PlotErrorGraph(res{3},'error_wg','d--','b','BarkerHendersonCutoff_2D, vap');
    
    PlotErrorGraph(res{4},'error_wl','>-','b','BarkerHendersonHardCutoff_2D, liq');
    PlotErrorGraph(res{4},'error_wg','>--','b','BarkerHendersonHardCutoff_2D, vap');
    
    PlotErrorGraph(res{5},'error_wl','<-','b','BarkerHenderson_2D, liq');
    PlotErrorGraph(res{5},'error_wg','<--','b','BarkerHenderson_2D, vap');
        
    PlotErrorGraph(res{6},'error_wl','^-','m','ExponentialDouble, liq');
    PlotErrorGraph(res{6},'error_wg','^--','m','ExponentialDouble, vap');
    
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);
    xlabel('$N$','Interpreter','Latex','fontsize',20);
    ylabel(['Relative sum rule error $\frac{n(0)-n_C}{n_C}$'],...
                            'Interpreter','Latex','fontsize',20);
    xlim([(N(1)-2),(N(end)+2)]);    
    %legend(legendstring,'Location','northeast');%,'Orientation','horizontal');
    legend(legendstring,'Location','eastOutside');%,'Orientation','horizontal');
    
    comment = 'Computed for hard wall, hard sphere fluid and BH fluid';    
    SaveFigure('SumRuleError',v2struct(N,NS,config,comment));
	    
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
                
                if(~isfield(conf.optsPhys,'V2'))                                        
                    [~,~,params] = CL.Compute1D(in.eta);
                    res(i,j).error_wl = params.contactDensity_relError;
                else
                    res(i,j).error_conv1 = preErr.error_conv1;
                    res(i,j).Conv        = CL.IntMatrV2.Conv;
                    
                    [~,~,params] = CL.Compute1D('WL');
                    res(i,j).error_wl = params.contactDensity_relError;

                    [~,~,params] = CL.Compute1D('WG');
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