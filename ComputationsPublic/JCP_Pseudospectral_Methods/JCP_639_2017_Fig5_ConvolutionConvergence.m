function JCP_639_2017_Fig5_ConvolutionConvergence()
    
    AddPaths('JCP_639_2017');

        
    V1       = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.865);
    V2       = struct('V2DV2','BarkerHendersonCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5); 
    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      
                  
    PhysArea  = struct('N',[20,20],...
                       'L1',5,'L2',4,'L2_AD',4.,...
                       'y2wall',0.,...
                       'N2bound',10,'h',1,...
                       'alpha_deg',90);                  
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',10,'N2disc',10);                   
    V2Num     = struct('Fex','SplitAnnulus','N',[10,10]);                    
    optsNum   = struct('PhysArea',PhysArea,...
                       'FexNum',Fex_Num,...
                       'V2Num',V2Num,...
                       'maxComp_y2',-1,...
                       'y1Shift',0);

    config = v2struct(optsNum,optsPhys);                                
                
    
    N    = 20; 
    NS_d = 10; 
    NS   = 10:NS_d:80; 
    
    res = DataStorage('ConvError',@ComputeError,v2struct(N,NS,config),[],[],{'config_optsNum_PhysArea_N','config_optsPhys_V1_epsilon_w'});
    
    
    %**************************
    %******** Plotting ********
    %**************************    
    figure('color','white','Position',[0 0 900 800]);     
    
    gray         = 0.5*[1 1 1];
    legendstring = {};    
    PlotMatrixError(res,'ConvSplitDisk','*','k'); legendstring(end+1) = {'\phi_{SD}'};        
    PlotMatrixError(res,'ConvBHCutoff','o','k');  legendstring(end+1) = {'\phi_{BH,rC}'};  
    
    PlotMatrixError(res,'A_n2','s','k'); legendstring(end+1) = {'w_2'}; 
    PlotMatrixError(res,'A_n3','p','k'); legendstring(end+1) = {'w_3'}; 
    PlotMatrixError(res,'A_n2_v_1','v','k'); legendstring(end+1) = {'w_{2,1}'}; 
    PlotMatrixError(res,'A_n2_v_2','<','k'); legendstring(end+1) = {'w_{2,2}'}; 
    
    PlotMatrixError(res,'AAD_n2','s',gray); legendstring(end+1) = {'w_{2,AD}'}; 
    PlotMatrixError(res,'AAD_n3','p',gray); legendstring(end+1) = {'w_{3,AD}'}; 
    PlotMatrixError(res,'AAD_n2_v_1','v',gray); legendstring(end+1) = {'w_{2AD,1}'}; 
    PlotMatrixError(res,'AAD_n2_v_2','<',gray); legendstring(end+1) = {'w_{2AD,2}'}; 
    
    legend(legendstring);
    set(gca,'YScale','log');
    set(gca,'linewidth',1.5);
    set(gca,'fontsize',20);
    xlabel('$M$','Interpreter','Latex','fontsize',20);
    ylabel(['$\|{\bf C}_{\phi}(M)-{\bf C}_\phi(M+',num2str(NS_d),')\|_2$'],...
                            'Interpreter','Latex','fontsize',20);
    xlim([(NS(1)-2),(NS(end-1)+2)]);
    ylim([1e-16 1]);
    SaveFigure('JCP_639_2017_Fig5',v2struct(N,NS,config));

    
    function res = ComputeError(in,h)
        conf   = in.config;        
        ns     = in.NS;        
                
        for i = 1:length(in.N)
            
            conf.optsNum.PhysArea.N = [in.N(i),in.N(i)]; 
            
            CL = ContactLineHS(conf);
            CL.Preprocess();        
            
            
            for j = 1:length(ns)
                CL.optsNum.V2Num.N       = [ns(j),ns(j)];
                CL.optsNum.FexNum.N1disc = ns(j);
                CL.optsNum.FexNum.N2disc = ns(j);                                
                
                res(i,j).N  = in.N(i);
                res(i,j).NS = ns(j); 
                
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

                CL.optsNum.V2Num.Fex            = 'SplitDisk';
                CL.optsPhys.V2.V2DV2            = 'BarkerHenderson_2D'; 
                preErr                          = CL.Preprocess_MeanfieldContribution();                
                res(i,j).error_conv_SplitDisk   = preErr.error_conv1;
                res(i,j).ConvSplitDisk          = CL.IntMatrV2.Conv;
        
                close all;                       
            end
            clear('CL');         
        end
    end   

    function PlotMatrixError(res,A_name,sym,col)
        
        n = 1;
        
        for k1 = 1:size(res,1)
            for k2 = 1:(size(res,2)-1)
                line(n)   = norm(res(k1,k2).(A_name)-res(k1,k2+1).(A_name));                
                line_N(n) = res(k1,k2).NS;
                n = n+1;
            end
        end     
        plot(line_N,line,['-',sym],'color',col,'MarkerSize',10,...
                                   'MarkerFaceColor',col); hold on;                   
    end

end