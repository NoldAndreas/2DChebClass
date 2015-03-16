function ComputeSurfaceTensions

    AddPaths('CodePaper');            
    
    PhysArea = struct('N',[20,20],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',90);
                      
    V2Num    = struct('Fex','SplitDisk','N',[80,80]);
    V2       = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',5);     

    FexNum   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);

    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',FexNum,'V2Num',V2Num,...
                     'maxComp_y2',20,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.0);%1.375);%1.25)s;
    

    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',0.75,...                                               
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);

                  
    config = v2struct(optsNum,optsPhys);      
    
    res = DataStorage('SurfaceTensionComputations',@ComputeSurfaceTensions,config,[]);
    
    function res = ComputeSurfaceTensions(config,h)
        %(1) Compute Liquid-Gas Surface Tensions
        config.optsNum.PhysArea.N = [80,3];
        config.optsNum.PhysArea.N2bound = 3;

        CL = ContactLineHS(config);
        CL.Preprocess();    

        CL.optsNum.V2Num.Fex   = 'SplitAnnulus';
        CL.optsPhys.V2.V2DV2   = 'BarkerHendersonCutoff_2D';                            
        preErr                 = CL.Preprocess_MeanfieldContribution();                
        res{1}.gammaLG         = CL.Compute1D('LG');
        res{1}.config          = CL.GetConfig();

        CL.optsNum.V2Num.Fex     = 'SplitDisk';
        CL.optsPhys.V2.V2DV2     = 'BarkerHenderson_2D'; 
        preErr         = CL.Preprocess_MeanfieldContribution();                
        res{2}.gammaLG = CL.Compute1D('LG');
        res{2}.config  = CL.GetConfig();

        CL.optsNum.V2Num.Fex     = 'SplitDisk';                
        CL.optsPhys.V2.V2DV2     = 'ExponentialDouble';
        preErr             = CL.Preprocess_MeanfieldContribution();                
        res{3}.gammaLG = CL.Compute1D('LG');
        res{3}.config  = CL.GetConfig();



        %(2) Compute Wall-Fluid Surface Tensions
        config.optsNum.PhysArea.N = [1,80];
        config.optsNum.PhysArea.N2bound = 20;

        CL = ContactLineHS(config);
        CL.Preprocess();    

        CL.optsNum.V2Num.Fex   = 'SplitAnnulus';
        CL.optsPhys.V2.V2DV2   = 'BarkerHendersonCutoff_2D';                            
        preErr                 = CL.Preprocess_MeanfieldContribution();                
        res{1}.gammaWL         = CL.Compute1D('WL');
        res{1}.gammaWG         = CL.Compute1D('WG');
        res{1}.config          = CL.GetConfig();

        CL.optsNum.V2Num.Fex   = 'SplitDisk';
        CL.optsPhys.V2.V2DV2   = 'BarkerHenderson_2D'; 
        preErr                 = CL.Preprocess_MeanfieldContribution();                
        res{2}.gammaWL         = CL.Compute1D('WL');
        res{2}.gammaWG         = CL.Compute1D('WG');    
        res{2}.config          = CL.GetConfig();

        CL.optsNum.V2Num.Fex     = 'SplitDisk';                
        CL.optsPhys.V2.V2DV2     = 'ExponentialDouble';
        preErr                   = CL.Preprocess_MeanfieldContribution();                
        res{3}.gammaWL           = CL.Compute1D('WL');
        res{3}.gammaWG           = CL.Compute1D('WG');    
        res{3}.config            = CL.GetConfig();

    end
end