function ThesisNanoscale_Fig3_ComputeContactAngles()

   global dirData
   AddPaths('ThesisNanoscale');            
            
   opts.bounds1   = [0,20];
   opts.alpha_deg = 45;  
   opts.epw       = 1.155;   
   opts.dryingWetting = 'wetting';   
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
   Job_ComputeContactAngle(opts);    
    
   opts.bounds1   = [-10,10];
   opts.alpha_deg = 60;  
   opts.epw       = 1.071;
   opts.dryingWetting = 'wetting';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);   
   Job_ComputeContactAngle(opts);                
   
   opts.alpha_deg = 90;  
   opts.epw       = 0.856;     
   opts.dryingWetting = 'wetting';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
   Job_ComputeContactAngle(opts);       
   
   opts.alpha_deg = 120;  
   opts.epw       = 0.594;
   opts.dryingWetting = 'drying';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
   Job_ComputeContactAngle(opts); 
   
   opts.alpha_deg = 135;  
   opts.epw       = 0.453;
   opts.dryingWetting = 'drying';
   opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
   Job_ComputeContactAngle(opts);

    function Job_ComputeContactAngle(opts)

        config = GetStandardConfig(opts);
        close all;

        CLT = ContactLineHS(config);     
        CLT.Preprocess();
        CLT.ComputeEquilibrium();      
        CLT.PostProcess(opts);
        CLT.PlotDensitySlices();
        CLT.PlotDisjoiningPressures();                
    end
    function config = GetStandardConfig(opts)
        
        alpha_deg  = opts.alpha_deg;
        epw        = opts.epw;
        
        bounds1    = opts.bounds1;
        bounds2    = [0.5,15.5];
        maxComp_y2 = 35;        
        N          = [50,80];           

        PhysArea = struct('N',N,'L1',4,'L2',2,...
                          'alpha_deg',alpha_deg);

        PlotArea = struct('y1Min',bounds1(1),'y1Max',bounds1(2),...
                          'y2Min',bounds2(1),'y2Max',bounds2(2),...
                          'zMax',4,...
                          'N1',100,'N2',100);

        V2Num     = struct('Fex','SplitAnnulus','N',[80,80]);          
        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid','Ncircle',1,'N1disc',50,'N2disc',50); 

        optsNum = struct('PhysArea',PhysArea,...
                         'PlotAreaCart',PlotArea,...
                         'FexNum',Fex_Num,...
                         'V2Num',V2Num,...
                         'maxComp_y2',maxComp_y2);

        V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',epw);
        V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5); 

        optsPhys = struct('V1',V1,'V2',V2,...                   
                          'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

        config = v2struct(optsNum,optsPhys);                        

    end
    function filename = ComputeExactAdsorptionIsotherm(opts)
       
        PhysArea = struct('N',[1,250],...
                          'L1',4,'L2',2,...                          
                          'alpha_deg',90);

        V2Num     = struct('Fex','SplitAnnulus','N',[80,80]);                  
        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid','Ncircle',1,'N1disc',50,'N2disc',50);        

        optsNum = struct('PhysArea',PhysArea,...
                         'FexNum',Fex_Num,...
                         'V2Num',V2Num,...
                         'maxComp_y2',-1,...
                         'y1Shift',0);

        V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',opts.epw);
        V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5); 

        optsPhys = struct('V1',V1,'V2',V2,...                   
                          'kBT',0.75,...                                                    
                          'Dmu',0.0,'nSpecies',1,...
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