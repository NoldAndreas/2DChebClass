function MMNP_Fig3_ComputeContactAngles()

    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'MMNP'],'ORG');    
            
   opts.bounds1   = [0,20];
   opts.alpha_deg = 40;  
   opts.epw       = 1.375;   
   opts.dryingWetting = 'wetting';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
   Job_ComputeContactAngle(opts);    
   
   opts.bounds1   = [-10,10];
   opts.alpha_deg = 60;  
   opts.epw       = 1.25;   
   opts.dryingWetting = 'wetting';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);   
   Job_ComputeContactAngle(opts);                
   
   opts.alpha_deg = 90;  opts.epw       = 1.0;     
   opts.dryingWetting = 'wetting';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
   Job_ComputeContactAngle(opts);       
   
   opts.alpha_deg = 120;  opts.epw      = 0.7;   
   opts.dryingWetting = 'drying';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
   Job_ComputeContactAngle(opts); 
   
   opts.alpha_deg = 135;  opts.epw      = 0.55;
   opts.dryingWetting = 'drying';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
   Job_ComputeContactAngle(opts);

    function Job_ComputeContactAngle(opts)

        config = GetStandardConfig(opts);
        close all;

        CLT = ContactLineHS(config);     
        CLT.Preprocess();
        CLT.ComputeEquilibrium();      
        CLT.PostProcess(opts.AdsorptionIsotherm_file);
        CLT.PlotDensitySlices();
        
        %CLT.PlotDensityResult();
        %CLT.PlotContourResults();
        %CLT.PlotDisjoiningPressures();

        %     CLT.ComputeAdsorptionIsotherm('load'); %load \2DChebData\POF_FMT_ContactLine\deg90\IterativeContinuationPostProcess\2014_1_20_18_46
    %     CLT.PostProcess_2DDisjoiningPressure();
    %     
    % 	%f1 = figure('Color','white','Position',[0 0 1000 600]);
    %     [f1,f2] = CLT.Post_HFrom2DDisjoiningPressure();
    %
    %     print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_Interfaces'],f1);
    %     saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_Interfaces.fig']);
    %
    %     print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_DisjoiningPressure'],f2);
    %     saveas(f2,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_DisjoiningPressure.fig']);
    %     
    %     inset2(f1,f2,0.4,[0.26,0.55]);
    %     %inset2(f1,f2,0.35,[0.22,0.55]);
    %     close(f2);      
    %     
    %     print2eps([dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure'],f1);
    %     saveas(f1,[dirData filesep 'EquilibriumSolutions' filesep CLT.FilenameEq '_InterfacesAndDisjoiningPressure.fig']);
    %         
    %     CLT.FittingAdsorptionIsotherm([10 14],1);
    %     CLT.SumRule_DisjoiningPotential();
    %     %***************************************  
    end

    function config = GetStandardConfig(opts)
        
        alpha_deg  = opts.alpha_deg;
        epw        = opts.epw;
        
        bounds1    = opts.bounds1;
        bounds2    = [0.5,15];
        maxComp_y2 = 35;        
        N          = [50,80];           

        PhysArea = struct('N',N,'L1',4,'L2',2,'L2_AD',2.,...
                          'y2wall',0.,...
                          'N2bound',14,'h',1,...
                          'alpha_deg',alpha_deg);

        PlotArea = struct('y1Min',bounds1(1),'y1Max',bounds1(2),...
                          'y2Min',bounds2(1),'y2Max',bounds2(2),...
                          'zMax',4,...
                          'N1',100,'N2',100);

        V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);    
        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                           'Ncircle',1,'N1disc',50,'N2disc',50); %35,34

        optsNum = struct('PhysArea',PhysArea,...
                         'PlotAreaCart',PlotArea,...
                         'FexNum',Fex_Num,...
                         'V2Num',V2Num,...
                         'maxComp_y2',maxComp_y2);

        V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',epw);
        V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

        optsPhys = struct('V1',V1,'V2',V2,...                   
                          'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

        config = v2struct(optsNum,optsPhys);                        

    end

    function filename = ComputeExactAdsorptionIsotherm(opts)
       
        PhysArea = struct('N',[1,250],...
                          'L1',5,'L2',4,'L2_AD',4.,...
                          'y2wall',0.,...
                          'N2bound',24,'h',1,...
                          'alpha_deg',90);

        V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);

        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                           'Ncircle',1,'N1disc',50,'N2disc',50);                   

        optsNum = struct('PhysArea',PhysArea,...
                         'FexNum',Fex_Num,...
                         'V2Num',V2Num,...
                         'maxComp_y2',-1,...
                         'y1Shift',0);

        V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',opts.epw);
        V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

        optsPhys = struct('V1',V1,'V2',V2,...                   
                          'kBT',0.75,...                                                    
                          'Dmu',0.03,'nSpecies',1,...
                          'sigmaS',1);      

        config = v2struct(optsNum,optsPhys);                        

        CL = ContactLineHS(config);
        CL.Preprocess();    close all;
        CL.ComputeAdsorptionIsotherm(1000,opts.dryingWetting);    
        
        CL.FittingAdsorptionIsotherm([10 14],1)
        if(optsPhys.kBT == 0.75)
            CL.SumRule_AdsorptionIsotherm(0.3463);
        end
        
        filename = CL.AdsorptionIsotherm.Filename;


    end

    
end