function CheckSumRule_HalfSpace_L(N)
    
    
    AddPaths('CodePaper');    
    
    if(nargin == 0)
        N = 100;
    end    

    PhysArea = struct('N',[1,N],...
                      'L1',5,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',16,'h',1,...
                      'alpha_deg',90);
                      

    V2Num   = struct('Fex','SplitDisk','N',[80,80]);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 
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

    L = 0.25:0.25:4;    
    folder = 'SumRuleErrorL';
        
    eta = 0.3;    
    res{1} = DataStorage(folder,@ComputeError,v2struct(L,config,eta),[]);
        	    
    eta = 0.15;    
    res{2} = DataStorage(folder,@ComputeError,v2struct(L,config,eta),[]);  
        
    config.optsNum.V2Num = struct('Fex','SplitAnnulus','N',[80,80]);
    config.optsPhys.V2   = struct('V2DV2','BarkerHendersonCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',5);
    config.optsPhys.V1.V1DV1 = 'Vext_BarkerHenderson_HardWall';
    config.optsPhys.V1.epsilon_w = 0.9; %0.95    
    res{3} = DataStorage(folder,@ComputeError,v2struct(L,config),[]);     	
   
%    config.optsNum.V2Num  = struct('Fex','SplitDisk','N',[80,80]); 
%    config.optsPhys.V2    = struct('V2DV2','ExponentialDouble','epsilon',1,'LJsigma',1);
%    config.optsPhys.V1.epsilon_w = 1.45;    
%    res{4}  = DataStorage(folder,@ComputeError,v2struct(L,config),[]);    
        
                    
    figure('color','white','Position',[0 0 900 800]); 
    legendstring = {};
	PlotErrorGraph(res{1},'error_wl','s-','k',['Hard sphere, eta = ',num2str(0.3)]); 
    PlotErrorGraph(res{2},'error_wl','p-','k',['Hard sphere, eta = ',num2str(0.15)]);     
    
    PlotErrorGraph(res{3},'error_wl','o-','k','BarkerHendersonCutoff_2D, liq');
    PlotErrorGraph(res{3},'error_wg','o--','k','BarkerHendersonCutoff_2D, vap');
        
%    PlotErrorGraph(res{4},'error_wl','^-','k','ExponentialDouble, liq');
%    PlotErrorGraph(res{4},'error_wg','^--','k','ExponentialDouble, vap');   
       
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);
    xlabel('$L$','Interpreter','Latex','fontsize',20);
    ylabel(['Relative sum rule error $\frac{n(0)-n_C}{n_C}$'],...
                            'Interpreter','Latex','fontsize',20);
    xlim([0,(L(end)+0.5)]);    
    
     set(gca,'YScale','log');
    %legend(legendstring,'Location','northeast');%,'Orientation','horizontal');
   % legend(legendstring,'Location','eastOutside');%,'Orientation','horizontal');
    
    comment = 'Computed for hard wall, hard sphere fluid and BH fluid';    
    SaveFigure('SumRuleError_L',v2struct(L,config,comment));
	    
   
    function res = ComputeError(in,h)
        conf = in.config;
        l      = in.L;        
        
        %error_wg    = zeros(length(N),length(NS));
        %error_wl    = zeros(length(N),length(NS));
                
        for i = 1:length(l)
            
            conf.optsNum.PhysArea.L2    = l(i);
            conf.optsNum.PhysArea.L2_AD = l(i);            
            
            CL = ContactLineHS(conf);
            preErr = CL.Preprocess(); 

            if(~isfield(conf.optsPhys,'V2'))                                        
                [~,~,params] = CL.Compute1D(in.eta);
                res(i).error_wl = params.contactDensity_relError;
            else
                res(i).error_conv1 = preErr.error_conv1;
                res(i).Conv        = CL.IntMatrV2.Conv;

                [~,~,params] = CL.Compute1D('WL');
                res(i).error_wl = params.contactDensity_relError;

                [~,~,params] = CL.Compute1D('WG');
                res(i).error_wg = params.contactDensity_relError;
            end
                
            res(i).L  = l(i);
            close all;
            clear('CL');                
        end
%        res.error_wl = error_wl;
%        res.error_wg = error_wg; 
    end     
    function PlotErrorGraph(res,var_name,sym,col,name)
        
        n = 1;
                
        for k1 = 1:length(res)
            line(n)   = abs(res(k1).(var_name));
            line_N(n) = (res(k1).L);%+res(k1,k2+1).NS)/2;
            %plot(line_N(n),line(n),...
             %       [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;
            n = n+1;
        end
        plot(line_N,line,...
                        [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;        
                    
        legendstring(end+1) = {name};
    end

end