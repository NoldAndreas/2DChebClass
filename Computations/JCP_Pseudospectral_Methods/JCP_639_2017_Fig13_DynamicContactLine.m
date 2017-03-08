function JCP_639_2017_Fig13_DynamicContactLine()

    AddPaths('JCP_639_2017');            
    close all;
    
	%**************************
    % **** Initialization *****
    %**************************
    
    PhysArea  = struct('N',[20,20],...
                       'L1',4,'L2',2,...                                            
                       'alpha_deg',90);
	SubArea   = struct('shape','Box','y1Min',-2,'y1Max',2,...
                          'y2Min',0.5,'y2Max',2.5,...
                          'N',[40,40]);                      
    V2Num     = struct('Fex','SplitAnnulus','N',[80,80]);
    FexNum    = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
	plotTimes = struct('t_int',[0,20],'t_n',100);
    optsNum   = struct('PhysArea',PhysArea,...
                       'FexNum',FexNum,'V2Num',V2Num,...
                       'SubArea',SubArea,...
                       'maxComp_y2',20,...
                       'y1Shift',0,...
                       'plotTimes',plotTimes);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.856,...
                'tau',5,'epsilon_w_max',1);
	V2       = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);     
    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);

    config = v2struct(optsNum,optsPhys);                                    
    N       = 20:10:70;    
    
    %************************
    % **** Computations *****
    %************************
    
    
    ignoreList = {'config_optsNum_PhysArea_N','config_optsNum_PlotAreaCart'};        
    res{1} = DataStorage('DynamicError',@ComputeDynamicError,v2struct(config,N),[],[],ignoreList);    
                
    config.optsPhys.Inertial = true;    
    config.optsPhys.gammaS   = 2;        
    res{2} = DataStorage('DynamicError',@ComputeDynamicError,v2struct(config,N),[],[],ignoreList);    
    
    res = PostProcess(res);    
   
	%********************
    % **** Plotting *****
    %********************   	                
    PlotResultsOverTime(res,'mass',{'maxN'}); ylim([1.9 2.6]); f1 = gcf; 
    PlotResults(res,'massErrorRel_dtMax');    f2 = gcf;  
    inset2(f2,f1,0.3,[0.27,0.2]);  close(f1);    
    SaveFigure('JCP_639_2017_Fig13_a');   
    
    PlotExampleSnaptshots(res);
    SaveFigure('JCP_639_2017_Fig13_b');
        
    function PlotExampleSnaptshots(res)
        conf = res{1}(3).config;
        conf.optsNum.PlotAreaCart = struct('y1Min',-5,'y1Max',5,...
                                           'y2Min',0.5,'y2Max',10.5,...
                                           'N1',100,'N2',100,...
                                           'NFlux',15);

        CL = ContactLineHS(conf);
        CL.Preprocess(); 
        CL.ComputeEquilibrium();              
        CL.ComputeDynamics();            
        CL.PostprocessDynamics();
            
        plotData.subArea  = CL.subArea;
        plotData.optsPhys = CL.optsPhys;
        plotData.optsNum  = CL.optsNum;
        plotData.data     = CL.dynamicsResult;
        plotData.filename = CL.FilenameDyn;
        plotData.data.shape = CL.IDC;
        plotData.data.t = plotData.data.t(plotData.data.t < 3.1);
        
        figure('color','white','Position',[0 0 800 800]);
        PlotDDFT_SnapshotsShape(plotData,[],{'4Snapshots','noNewFigure','NumericsManuscript'});                        
    end
    function res = PostProcess(res)
                
        for i0 = 1:length(res)            
            NT = length(res{i0}(1).t);
            T  = res{i0}(1).t(end);
            [xt,w,D] = FDxw(NT);
            if(max(abs((xt+1)*T/2 - res{i0}(1).t')) < 1e-12)
                Dt = D*2/T;
            else
                disp('barychebdiff used');
                Diff_T = barychebdiff(res{i0}(1).t',1);
                Dt = Diff_T.Dx;
            end
            
            for i1 = 1:length(res{i0})
                res{i0}(i1).massError_dt = Dt*res{i0}(i1).massError';
                res{i0}(i1).massError_dtMax = max(abs(res{i0}(i1).massError_dt));
                
                res{i0}(i1).massErrorRel_dt = Dt*res{i0}(i1).massErrorRel';
                res{i0}(i1).massErrorRel_dtMax = max(abs(res{i0}(i1).massErrorRel_dt));
            end
        end
    end    
    function res = ComputeDynamicError(opts,h)

        conf = opts.config;
        n      = opts.N;
        
        for i = 1:length(n)
            conf.optsNum.PhysArea.N       = n(i)*[1,1];
            
            conf.optsNum.FexNum.N1disc    =  2*round((20+n(i)/4)/2)*[1,1];
            conf.optsNum.FexNum.N2disc    =  2*round((20+n(i)/4)/2)*[1,1];
            conf.optsNum.V2Num.N          =  2*round((20+n(i)/2)/2)*[1,1];
            
            conf.optsNum.SubArea.N        =  n(i)*[1,1];

            CL = ContactLineHS(conf);
            CL.Preprocess(); 
            CL.ComputeEquilibrium(struct('solver','Picard'));              
            CL.ComputeDynamics();            
            CL.PostprocessDynamics();
       
            res(i).config          = conf;
            res(i).t               = CL.dynamicsResult.t;
            res(i).mass            = CL.dynamicsResult.Subspace.mass;
            res(i).massError       = CL.dynamicsResult.Subspace.massError;
            res(i).massErrorRel    = CL.dynamicsResult.Subspace.massErrorRel;
            res(i).MaxmassError    = max(abs(CL.dynamicsResult.Subspace.massError));
            res(i).MaxmassErrorRel = max(abs(CL.dynamicsResult.Subspace.massErrorRel));

            close all;
            clear('CL');            
        end        
    end
    function PlotResultsOverTime(res,var,opts)
        if(nargin < 3)
            opts = {};
        end
        lines = {'-','--',':','-.'};     
        cols = {}; 
        nocols = length(res{1});
        for iC = 1:nocols
            cols{end+1} = (nocols-iC)/nocols*[1 1 1];
        end 
        
        
        figure('color','white','Position',[0 0 800 800]); 
        
        for i0 = 1:length(res)
            lin = lines{i0};
            if(IsOption(opts,'maxN'))
                plot(res{i0}(end).t,res{i0}(end).(var),[lin],'linewidth',1.5,'color','k'); hold on;                
            else
                for i = 1:length(res{i0})            
                    plot(res{i0}(i).t,res{i0}(i).(var),[lin],'linewidth',1.5,'color',cols{i}); hold on;
                    n(i)  = res{i0}(i).config.optsNum.PhysArea.N(1);
                end
            end
        end
                
        ylim([2.2 3]);
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);            
        xlabel('$t$','Interpreter','Latex','fontsize',20);
        ylabel(GetYLabel(var),'Interpreter','Latex','fontsize',20);
                                       
    end
    function PlotResults(res,var)
        lines = {'-','--',':','-.'};     
        figure('color','white','Position',[0 0 800 800]);                 
        
        for i0 = 1:length(res)
            
            lin = lines{i0};
            n   = zeros(1,length(res{i0}));
            err = zeros(1,length(res{i0}));
            
            for i = 1:length(res{i0})
                n(i)   = res{i0}(i).config.optsNum.PhysArea.N(1);
                err(i) = res{i0}(i).(var);
            end
            
            plot(n,err,['ko',lin],'linewidth',1.5); hold on;
        end        
        set(gca,'YScale','log');
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',20);            
        xlabel('$N$','Interpreter','Latex','fontsize',20);
        ylabel(GetYLabel(var),...
                                'Interpreter','Latex','fontsize',20);       
    end
    function str = GetYLabel(var)        

        if(strcmp(var,'mass'))
            str = '$m$';
        elseif(strcmp(var,'massErrorRel_dt'))
            str = '$\frac{d}{dt} \left( \frac{m_{err}}{m} \right)$';
        elseif(strcmp(var,'massErrorRel_dtMax'))
            str = '$\max |\frac{d}{dt} \left( \frac{m_{err}}{m} \right)|$';
        elseif(strcmp(var,'massError_dtMax'))
            str = '$\max|\frac{d m_{err}}{dt}|$';
        elseif(strcmp(var,'massError_dt'))
            str = '$\frac{d m_{err}}{dt}$';
        elseif(strcmp(var,'MaxmassErrorRel'))
            str = 'Relative error of Mass in Subspace';
        elseif(strcmp(var,'MaxmassError'))
            str = 'Error of Mass in Subspace';
        elseif(strcmp(var,'massErrorRel'))
            str = 'Relative error of Mass in Subspace';
        elseif(strcmp(var,'massError'))       
            str = 'Error Mass in Subspace';
        end
    end
end