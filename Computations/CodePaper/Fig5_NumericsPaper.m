function CheckConvolutionHalfSpace_BH_Conv1()
    
    AddPaths('CodePaper');

    PhysArea = struct('N',[20,20],...
                      'L1',5,'L2',4,'L2_AD',4.,...
                      'y2wall',0.,...
                      'N2bound',10,'h',1,...
                      'alpha_deg',90);

    %V2Num   = struct('Fex','SplitDisk','N',[30,30]);
    V2Num   = struct('Fex','SplitAnnulus','N',[10,10]);
    V2 = struct('V2DV2','BarkerHendersonCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5); %r_cutoff = 5
%    V2Num   = struct('Fex','ConstShortRange','N',[30,30]);


    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',10,'N2disc',10);
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'V2Num',V2Num,...
                     'maxComp_y2',-1,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.7);%1.375);%1.25)s;
    %V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1);     
%    V2 = struct('V2DV2','ConstShortRange','epsilon',1,'LJsigma',1,'lambda',1.5); 
%    V2 = struct('V2DV2','Exponential','epsilon',1.5,'LJsigma',1); 
%    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1.,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %epw = 0.9;%[0.75,0.8,0.85,0.9,0.95];
    config.optsPhys.V1.epsilon_w = 0.865;%0.9;%    1.0;%1.25;%0.55;% 1.375; %0.7;%1.25;%375;%25; %375;%47;%1.25;
                
    N    = 20;%:10:50;    
    NS_d = 10; %2;  %10;
    NS   = 10:NS_d:80; %[20,30,40,50,60,70,80,82];%10:NS_d:50;%10:10:40;[20,22,30,32];
    
    res = DataStorage(['ConvError'],@ComputeError,v2struct(N,NS,config),[],[],{'config_optsNum_PhysArea_N','config_optsPhys_V1_epsilon_w'});
    
	%CLT = ContactLineHS(config);
	%CLT.PreprocessIDC();     
    
    nocols = 5;
    nosyms = 5;
    cols = {'b','m','k','r','g'};
    syms = {'d','s','o','>','<'};
    
    %**********************************************        
    %PlotMatrixErrorOverY_OnePoint('Conv',1);
    %PlotMatrixErrorOverY_OnePoint('Conv',15);
	%PlotMatrixErrorOverY('Conv');%res(j1,j2).A_n2-res(j1,j2+1).A_n2,syms{j2},cols{j2},['NS = ',num2str(res(j1,j2).NS)]);            
    %xlim([0,10]);
    
%    PlotMatrixErrorOverY('AAD_n3');%res(j1,j2).A_n2-res(j1,j2+1).A_n2,syms{j2},cols{j2},['NS = ',num2str(res(j1,j2).NS)]);            
%    xlim([0,10]);    
    %**********************************************
    gray = 0.5*[1 1 1];
    figure('color','white','Position',[0 0 900 800]);     
    legendstring = {};    
    %PlotMatrixError('ConvBH_HardCutoff','>','k','\phi_{HC}');
    %PlotMatrixError('ConvSplitDiskExp','^','k','\phi_{Exp}');
    PlotMatrixError('ConvSplitDisk','*','k','\phi_{SD}');
    %PlotMatrixError('ConvShortRange','d','k','\phi_{SR}');
    %PlotMatrixError('Conv','o','k','\phi_{LR}');
    PlotMatrixError('ConvBHCutoff','o','k','\phi_{BH,rC}');    
    
    PlotMatrixError('A_n2','s','k','w_2'); %    PlotMatrixError('A_n2','d','m','w_2'); 
    PlotMatrixError('A_n3','p','k','w_3');
    PlotMatrixError('A_n2_v_1','v','k','w_{2,1}'); %PlotMatrixError('A_n2_v_1','<','b','w_{2,1}');
    PlotMatrixError('A_n2_v_2','<','k','w_{2,2}'); %PlotMatrixError('A_n2_v_2','>','r','w_{2,2}');
    
    PlotMatrixError('AAD_n2','s',gray,'w_{2,AD}'); %PlotMatrixError('AAD_n2','*','m','w_{2,AD}');
    PlotMatrixError('AAD_n3','p',gray,'w_{3,AD}'); %PlotMatrixError('AAD_n3','x','g','w_{3,AD}');
    PlotMatrixError('AAD_n2_v_1','v',gray,'w_{2AD,1}');
    PlotMatrixError('AAD_n2_v_2','<',gray,'w_{2AD,2}');%PlotMatrixError('AAD_n2_v_2','v','r','w_{2AD,2}');
    
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);
    xlabel('$M$','Interpreter','Latex','fontsize',20);
    ylabel(['$\|{\bf C}_{\phi}(M)-{\bf C}_\phi(M+',num2str(NS_d),')\|_2$'],...
                            'Interpreter','Latex','fontsize',20);
    xlim([(NS(1)-2),(NS(end-1)+2)]);
    ylim([1e-16 1]);
    SaveFigure('ConvolutionError1',v2struct(N,NS,config));
    %legend(legendstring,'Location','northOutside','Orientation','horizontal');
   
%     filename = 'ConvolutionMatrixError';
% 	print2eps([dirData filesep filename],gcf);
% 	saveas(gcf,[dirData filesep filename '.fig']);        
%     disp(['Figures saved in ',dirData filesep filename '.fig/eps']);    
    
    figure('color','white','Position',[0 0 800 800]); 
%    PlotErrorGraph('error_convBH_HardCutoff','>','k');
%    PlotErrorGraph('error_conv_SplitDiskExp','^','k');
    PlotErrorGraph('error_conv_SplitDisk','*','k');
    %PlotErrorGraph('error_convShortRange','d','k');
%    PlotErrorGraph('error_conv1','o','k');

    PlotErrorGraph('error_n2_1','s','g');
    PlotErrorGraph('error_n3_1','p','g');
    PlotErrorGraph('error_n2v2_1','v','g');        
    PlotErrorGraph('error_n2AD_ones','s','b');
    PlotErrorGraph('error_n3AD_ones','p','b');
    PlotErrorGraph('error_n2v2AD_ones','v','b');    
    
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',15);
    xlabel('$N$','Interpreter','Latex','fontsize',15);
    ylabel('error$(\Phi_{2D} \ast {\bf 1})$','Interpreter','Latex','fontsize',15);
    xlim([(NS(1)-2),(NS(end-1)+2)]);
    ylim([1e-16 1]);
    
    SaveFigure('ConvolutionError2',v2struct(N,NS,config));
    
    function res = ComputeError(in,h)
        conf = in.config;
        N      = in.N;
        NS     = in.NS;
        
        %error_wg    = zeros(length(N),length(NS));
        %error_wl    = zeros(length(N),length(NS));
                
        for i = 1:length(N)
            
            conf.optsNum.PhysArea.N = [N(i),N(i)]; %N(i),N(i)
            
            CL = ContactLineHS(conf);
            preErr = CL.Preprocess();        
            
            
            for j = 1:length(NS)
                CL.optsNum.V2Num.N       = [NS(j),NS(j)];
                CL.optsNum.FexNum.N1disc = NS(j);
                CL.optsNum.FexNum.N2disc = NS(j);                                
                
                res(i,j).N  = N(i);
                res(i,j).NS = NS(j); 
                
                preErr = CL.Preprocess_HardSphereContribution();
                               
                res(i,j).error_n2_1    = preErr.error_n2_1;
                res(i,j).error_n3_1    = preErr.error_n3_1;
                res(i,j).error_n2v2_1  = preErr.error_n2v2_1;                                
                                
                res(i,j).error_n2AD_erf   = preErr.error_n2AD_erf;
                res(i,j).error_n3AD_erf   = preErr.error_n3AD_erf;
                res(i,j).error_n2v2AD_erf = preErr.error_n2v2AD_erf;
                
                res(i,j).error_n2AD_ones   = preErr.error_n2AD_ones;
                res(i,j).error_n3AD_ones   = preErr.error_n3AD_ones;
                res(i,j).error_n2v2AD_ones = preErr.error_n2v2AD_ones;
                                
                
                res(i,j).A_n2        = CL.IntMatrFex.AD.n2;
                res(i,j).A_n3        = CL.IntMatrFex.AD.n3;
                res(i,j).A_n2_v_1    = CL.IntMatrFex.AD.n2_v_1;
                res(i,j).A_n2_v_2    = CL.IntMatrFex.AD.n2_v_2;                
                res(i,j).AAD_n2      = CL.IntMatrFex.AAD.n2;
                res(i,j).AAD_n3      = CL.IntMatrFex.AAD.n3;
                res(i,j).AAD_n2_v_1  = CL.IntMatrFex.AAD.n2_v_1;
                res(i,j).AAD_n2_v_2  = CL.IntMatrFex.AAD.n2_v_2;
                
                CL.optsNum.V2Num.Fex            = 'SplitAnnulus';
                CL.optsPhys.V2.V2DV2            = 'BarkerHenderson_2D'; 
                preErr                          = CL.Preprocess_MeanfieldContribution();                
                res(i,j).error_conv_BHcutoff    = preErr.error_conv1;
                res(i,j).ConvBHCutoff           = CL.IntMatrV2.Conv;

%                 CL.optsNum.V2Num.Fex   = 'SplitAnnulus';
%                 CL.optsPhys.V2.V2DV2   = 'BarkerHendersonCutoff_2D';                            
%                 preErr                 = CL.Preprocess_MeanfieldContribution();                
%                 res(i,j).error_conv1   = preErr.error_conv1; 
%                 res(i,j).Conv          = CL.IntMatrV2.Conv;                
%                 
%                 CL.optsNum.V2Num.Fex   = 'SplitAnnulus'; 
%                 CL.optsPhys.V2.V2DV2   = 'BarkerHendersonHardCutoff_2D';
%                 preErr                 = CL.Preprocess_MeanfieldContribution();                
%                 res(i,j).error_convBH_HardCutoff   = preErr.error_conv1;
%                 res(i,j).ConvBH_HardCutoff         = CL.IntMatrV2.Conv;
                
                
%                 CL.optsNum.V2Num.Fex   = 'ConstShortRange';
%                 CL.optsPhys.V2.V2DV2   = 'ConstShortRange';
%                 CL.optsPhys.V2.lambda  = 1.5;
%                 preErr                 = CL.Preprocess_MeanfieldContribution();                
%                 res(i,j).error_convShortRange   = preErr.error_conv1;
%                 res(i,j).ConvShortRange         = CL.IntMatrV2.Conv;
%                 
                
                CL.optsNum.V2Num.Fex            = 'SplitDisk';
                CL.optsPhys.V2.V2DV2            = 'BarkerHenderson_2D'; 
                preErr                          = CL.Preprocess_MeanfieldContribution();                
                res(i,j).error_conv_SplitDisk   = preErr.error_conv1;
                res(i,j).ConvSplitDisk          = CL.IntMatrV2.Conv;
%                 
%                 CL.optsNum.V2Num.Fex               = 'SplitDisk';                
%                 CL.optsPhys.V2.V2DV2               = 'ExponentialDouble';
%                 preErr                             = CL.Preprocess_MeanfieldContribution();                
%                 res(i,j).error_conv_SplitDiskExp   = preErr.error_conv1;
%                 res(i,j).ConvSplitDiskExp          = CL.IntMatrV2.Conv;
                
               % [~,~,params] = CL.Compute1D('WL');
               % res(i,j).error_wl = params.contactDensity_relError;

                %[~,~,params] = CL.Compute1D('WG');
                %res(i,j).error_wg = params.contactDensity_relError;
                
                close all;                       
            end
            clear('CL');         
        end
%        res.error_wl = error_wl;
%        res.error_wg = error_wg; 
    end   
    function PlotMatrixErrorOverY_OnePoint(A_name,kP)
        figure('color','white','Position',[0 0 800 800]);         
        leg_string = {};
        if(~isempty(strfind(A_name,'AAD')) || strcmp(A_name,'Conv'))
            y2 = CLT.IDC.Pts.y2_kv;
        else
            y2 = CLT.IDC.AD.Pts.y2_kv;
        end
        
        n = 1;
        for k1 = 1%:size(res,1)
            for k2 = 1:(size(res,2)-1)
                A = res(k1,k2).(A_name)-res(k1,k2+1).(A_name);
                plot(y2,abs(A(kP,:)),['-',syms{1+mod(n,nosyms)},...
                                         cols{1+mod(n,nocols)}],...
                                         'MarkerSize',10,...
                                         'MarkerFaceColor',cols{1+mod(n,nocols)}); hold on;
                leg_string(end+1) = {['NS = ',num2str(res(k1,k2).NS)]};
                n = n + 1;
            end
        end       
        xlim([0 10]);
        xlabel('$y_2$','Interpreter','Latex','fontsize',15);
        ylabel(['$|(',A_name,'_{N,y2_K}-',A_name,'_{N+',num2str(NS_d),',y2_K})|$'],...
                                'Interpreter','Latex','fontsize',15);
        title([A_name,'for y_2 = ',num2str(y2(kP))]);
        set(gca,'YScale','log');
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',15);
        legend(leg_string,'Location','southoutside','Orientation','horizontal');        
    end
    function PlotMatrixErrorOverY(A_name)
        figure('color','white','Position',[0 0 800 800]);         
        leg_string = {};
        if(~isempty(strfind(A_name,'AAD')) || strcmp(A_name,'Conv'))
            y2 = CLT.IDC.Pts.y2_kv;
        else
            y2 = CLT.IDC.AD.Pts.y2_kv;
        end
        
        n = 1;
        for k1 = 1%:size(res,1)
            for k2 = 1:(size(res,2)-1)
                A = res(k1,k2).(A_name)-res(k1,k2+1).(A_name);
                plot(y2,max(A,[],2),['-',syms{1+mod(n,nosyms)},...
                                         cols{1+mod(n,nocols)}],...
                                         'MarkerSize',10,...
                                         'MarkerFaceColor',cols{1+mod(n,nocols)}); hold on;
                leg_string(end+1) = {['NS = ',num2str(res(k1,k2).NS)]};
                n = n + 1;
            end
        end       
        xlabel('$y_2$','Interpreter','Latex','fontsize',15);
        ylabel(['$\max(',A_name,'_{N,y2}-',A_name,'_{N+',num2str(NS_d),',y2})$'],...
                                'Interpreter','Latex','fontsize',15);
        title(A_name);
        set(gca,'YScale','log');
        set(gca,'linewidth',1.5);
        set(gca,'fontsize',15);
        legend(leg_string,'Location','southoutside','Orientation','horizontal');        
    end
    function PlotMatrixError(A_name,sym,col,name)
        
        n = 1;
        
        for k1 = 1%:size(res,1)
            for k2 = 1:(size(res,2)-1)
                line(n)   = norm(res(k1,k2).(A_name)-res(k1,k2+1).(A_name));
                %line(n)   = max(max(abs(res(k1,k2).(A_name)-res(k1,k2+1).(A_name))));
                line_N(n) = (res(k1,k2).NS);%+res(k1,k2+1).NS)/2;                
                n = n+1;
            end
        end     
        plot(line_N,line,...
                        ['-',sym],'color',col,'MarkerSize',10,'MarkerFaceColor',col); hold on;
        %plot(line_N,line,col);        
        
        legendstring(end+1) = {name};
        
    end
    function PlotErrorGraph(var_name,sym,col)
        
        n = 1;
                
        for k1 = 1:size(res,1)
            for k2 = 1:size(res,2)
                line(n)   = res(k1,k2).(var_name);
                line_N(n) = (res(k1,k2).NS);%+res(k1,k2+1).NS)/2;
                %plot(line_N(n),line(n),...
                 %       [sym,col],'MarkerSize',10,'MarkerFaceColor',col); hold on;
                n = n+1;
                
                                
            end
        end
        plot(line_N,line,...
                        ['-',sym],'color',col,'MarkerSize',10,'MarkerFaceColor',col); hold on;        
    end

end