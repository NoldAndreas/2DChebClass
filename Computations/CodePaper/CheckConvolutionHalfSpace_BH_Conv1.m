function CheckConvolutionHalfSpace_BH_Conv1()

    global dirData

    PhysArea = struct('N',[1,60],...
                      'L1',5,'L2',4,'L2_AD',4.,...
                      'y2wall',0.,...
                      'N2bound',24,'h',1,...
                      'alpha_deg',90);

    V2Num   = struct('Fex','SplitDisk','N',[30,30]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'V2Num',V2Num,...
                     'maxComp_y2',-1,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.7);%1.375);%1.25)s;
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 
%    V2 = struct('V2DV2','Exponential','epsilon',1.5,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %epw = 0.9;%[0.75,0.8,0.85,0.9,0.95];
    config.optsPhys.V1.epsilon_w = 0.9;%    1.0;%1.25;%0.55;% 1.375; %0.7;%1.25;%375;%25; %375;%47;%1.25;
                
    N    = 10;%:10:50;
    NS_d = 8;
    NS   = 10:NS_d:50;%10:10:40;
    
    res = DataStorage([],@ComputeError,v2struct(N,NS,config),[]);
    
    function res = ComputeError(in,h)
        config = in.config;
        N      = in.N;
        NS     = in.NS;
        
        %error_wg    = zeros(length(N),length(NS));
        %error_wl    = zeros(length(N),length(NS));
                
        for i = 1:length(N)
            config.optsNum.PhysArea.N = [1,N(i)];
            for j = 1:length(NS)
                config.optsNum.V2Num.N       = [NS(j),NS(j)];
                config.optsNum.FexNum.N1disc = NS(j);
                config.optsNum.FexNum.N2disc = NS(j);

                CL = ContactLineHS(config);
                preErr = CL.Preprocess(); 
                
                res(i,j).error_conv1   = preErr.error_conv1;
                res(i,j).error_n2_1    = preErr.error_n2_1;
                res(i,j).error_n3_1    = preErr.error_n3_1;
                res(i,j).error_n2v2_1  = preErr.error_n2v2_1;                                
                
                res(i,j).N  = N(i);
                res(i,j).NS = NS(j); 
                %err         = CheckMeanfieldConvolution(CL);
                %%res(i,j).error_conv1 = err.error_conv1; 
                res(i,j).Conv        = CL.IntMatrV2.Conv;
                
                res(i,j).A_n2       = CL.IntMatrFex.AD.n2;
                res(i,j).A_n3       = CL.IntMatrFex.AD.n3;
                res(i,j).A_n2_v_1   = CL.IntMatrFex.AD.n2_v_1;
                res(i,j).A_n2_v_2   = CL.IntMatrFex.AD.n2_v_2;
                
                %[~,~,params] = CL.Compute1D('WL');
                %error_wl = params.contactDensity_relError;

                %[~,~,params] = CL.Compute1D('WG');
                %error_wg = params.contactDensity_relError;

                clear('CL');
                close all;
            end
        end
%        res.error_wl = error_wl;
%        res.error_wg = error_wg; 
    end     

    figure('color','white','Position',[0 0 800 800]); 
    PlotMatrixError('Conv','o','k');
    PlotMatrixError('A_n2','d','m');
    PlotMatrixError('A_n3','s','g');
    PlotMatrixError('A_n2_v_1','<','b');
    PlotMatrixError('A_n2_v_2','>','r');
    
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',15);
    xlabel('$N$','Interpreter','Latex','fontsize',15);
    ylabel(['$\|A_N-A_{N+',num2str(NS_d),'}\|$'],...
                            'Interpreter','Latex','fontsize',15);
    xlim([(NS(1)-2),(NS(end-1)+2)]);
    
   
    filename = 'ConvolutionMatrixError';
	print2eps([dirData filesep filename],gcf);
	saveas(gcf,[dirData filesep filename '.fig']);        
    disp(['Figures saved in ',dirData filesep filename '.fig/eps']);
    
    PlotErrorGraph('error_conv1');
    %Produce plots
    
    
%     figure('color','white','Position',[0 0 800 800]); 
%     for k1 = 1:size(res,1)
%         for k2 = 1:size(res,2)
%             plot(res(k1,k2).NS,res(k1,k2).errorData.error_conv1_posy2,'o','MarkerSize',7,'MarkerFaceColor','k'); hold on;
%         end
%     end
    
    function PlotMatrixError(A_name,sym,col)
        
        n = 1;
        
        for k1 = 1%:size(res,1)
            for k2 = 1:(size(res,2)-1)
                line(n)   = norm(res(k1,k2).(A_name)-res(k1,k2+1).(A_name));
                line_N(n) = (res(k1,k2).NS);%+res(k1,k2+1).NS)/2;
                plot(line_N(n),line(n),...
                        [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;
                n = n+1;
            end
        end        
        plot(line_N,line,col);        
        
    end

    function PlotErrorGraph(var_name)
        figure('color','white','Position',[0 0 800 800]); 
        for k1 = 1:size(res,1)
            for k2 = 1:size(res,2)
                plot(res(k1,k2).NS,res(k1,k2).(var_name),...
                    'o','MarkerSize',7,'MarkerFaceColor','k'); hold on;
            end
        end
        set(gca,'YScale','log');
    end

end