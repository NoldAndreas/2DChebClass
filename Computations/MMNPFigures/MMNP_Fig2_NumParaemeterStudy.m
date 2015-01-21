function MMNP_Fig2_NumParaemeterStudy()  
    
    
	alpha_deg = 135;
	epw       = 0.55;    
	bounds1   = [-15 10];
	bounds2   = [0.5 18];        
    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'MMNP'],'ORG');    
         

    PhysArea = struct('N',[45,75],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',alpha_deg);
                  
    PlotArea = struct('y1Min',bounds1(1),'y1Max',bounds1(2),...
                      'y2Min',bounds2(1),'y2Max',bounds2(2),...
                      'N1',100,'N2',100);                  

	V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);    
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50); %35,34
 
    optsNum = struct('PhysArea',PhysArea,...
                     'PlotAreaCart',PlotArea,...
                     'FexNum',Fex_Num,...
                     'V2Num',V2Num,...
                     'maxComp_y2',35);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',epw);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                                
    
    %Setup result file for this Job
    filename   = ([dirData filesep]); 
    filename   = [filename,'Job_MeasureContactAngles_epw_',...
                                getTimeStr(),'.txt'];
    
    Struct2File(filename,config,['Computed at ',datestr(now)]);    
    
    close all;    
	    
    CLT = ContactLineHS(config);
    CLT.Preprocess();
    %CLT.ComputeEquilibrium();    
    %**************************************        
    
    [y2_20,theta_20] = GetY2Theta(20);
    [y2_25,theta_25] = GetY2Theta(25);
    [y2_30,theta_30] = GetY2Theta(30);
    [y2_35,theta_35] = GetY2Theta(35);
    
    clear 'CLT'
    config2 = config;
    config2.optsNum.PhysArea.alpha_deg = 134.2;
    CLT2 = ContactLineHS(config2);
    CLT2.Preprocess();    
    [y2_20_134_2,theta_20_134_2] = GetY2Theta140(20);
    [y2_25_134_2,theta_25_134_2] = GetY2Theta140(25);
    [y2_30_134_2,theta_30_134_2] = GetY2Theta140(30);
    [y2_35_134_2,theta_35_134_2] = GetY2Theta140(35);
    
    clear 'CLT'
    clear 'CLT2'
    
    config2 = config;
    config2.optsNum.PhysArea.alpha_deg = 134;
    CLT2 = ContactLineHS(config2);
    CLT2.Preprocess();    
    [y2_20_140,theta_20_140] = GetY2Theta140(20);
    [y2_25_140,theta_25_140] = GetY2Theta140(25);
    [y2_30_140,theta_30_140] = GetY2Theta140(30);
    [y2_35_140,theta_35_140] = GetY2Theta140(35);
        
    %Load configuration \Config_2014_8_18_11_32_15.686.mat    
    f1 = figure('Color','white','Position',[0 0 900 800]);    
    
    plot(y2_20,theta_20,'m:','linewidth',1.5); hold on; 
    plot(y2_25,theta_25,'m-.','linewidth',1.5); hold on; 
    plot(y2_30,theta_30,'m--','linewidth',1.5); hold on; 
    plot(y2_35,theta_35,'m','linewidth',1.5); hold on; 
    
    plot(y2_20_134_2,theta_20_134_2,'k:','linewidth',1.5); hold on; 
    plot(y2_25_134_2,theta_25_134_2,'k-.','linewidth',1.5); hold on; 
    plot(y2_30_134_2,theta_30_134_2,'k--','linewidth',1.5); hold on; 
    plot(y2_35_134_2,theta_35_134_2,'k','linewidth',1.5); hold on; 
        
    plot(y2_20_140,theta_20_140,'b:','linewidth',1.5); hold on; 
    plot(y2_25_140,theta_25_140,'b-.','linewidth',1.5); hold on; 
    plot(y2_30_140,theta_30_140,'b--','linewidth',1.5); hold on; 
    plot(y2_35_140,theta_35_140,'b','linewidth',1.5); hold on; 
    
 %   [y2_90deg,theta_90deg] =  Load90DegComputation();
    
    %plot(y2_90deg,theta_90deg,'r','linewidth',1.5);
    
    xlim([5 37]);
    %ylim([60 60.8]);   
    %set(gca,'YTick',60:0.2:60.8);
    plot([5 40],CLT2.alpha_YCA*180/pi*[1 1],'r--','linewidth',2);
    plot([5 40],134*[1 1],'b--','linewidth',2);
    plot([5 40],134.2*[1 1],'k--','linewidth',2);
    plot([5 40],135*[1 1],'m--','linewidth',2);
    
    xlabel('$y/\sigma$','Interpreter','Latex','fontsize',25);
    ylabel('$\theta [^\circ]$','Interpreter','Latex','fontsize',25);
    set(gca,'linewidth',1.5,'fontsize',20);    

    f2 = figure('Color','white');
    CLT2.PlotContourResults(true);
        
    inset2(f1,f2,0.5,[0.45,0.55]);
    close(f2);   
    
    print2eps([dirData filesep 'CA_Asymptotics' filesep 'CA'],f1);
    saveas(f1,[dirData filesep 'CA_Asymptotics' filesep 'CA.fig']);
        
    %**************************************
     
%     configT.optsNum.maxComp_y2 = 25;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium();  
%     
%     configT.optsNum.maxComp_y2 = 30;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium(); 
%     
%     configT.optsNum.maxComp_y2 = 35;
%     CLT = ContactLine(configT);
%     CLT.Preprocess();
%     CLT.ComputeEquilibrium(); 
    
	function [y2,theta] = GetY2Theta(y2Max)
        CLT.optsNum.maxComp_y2 = y2Max;
        CLT.ComputeEquilibrium();  
        [y2,theta] = CLT.PlotInterfaceAnalysisY2([5 (y2Max+3)]);        
    end

    function [y2,theta] = GetY2Theta140(y2Max)
        CLT2.optsNum.maxComp_y2 = y2Max;
        CLT2.ComputeEquilibrium();  
        [y2,theta] = CLT2.PlotInterfaceAnalysisY2([5 (y2Max+3)]);        
    end

    function [y2,theta] = Load90DegComputation()

        PhysArea = struct('N',[50,80],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,'N2bound',14,'h',1,...
                      'alpha_deg',90);
    
        V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);    
        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                           'Ncircle',1,'N1disc',50,'N2disc',50);                   

        optsNum = struct('PhysArea',PhysArea,...
                         'FexNum',Fex_Num,...
                         'maxComp_y2',15,...
                         'V2Num',V2Num);

        V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.55);
        V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

        optsPhys = struct('V1',V1,'V2',V2,...                   
                          'kBT',0.75,'Dmu',0.0,...
                          'nSpecies',1,...
                          'sigmaS',1);      

        config = v2struct(optsNum,optsPhys);                        

        
        ChangeDirData();            
        
        CLT90 = ContactLineHS(config);     
        CLT90.Preprocess();
        CLT90.ComputeEquilibrium();       
        
        [y2,theta] = CLT90.PlotInterfaceAnalysisY2([5 20]);        
    end
end