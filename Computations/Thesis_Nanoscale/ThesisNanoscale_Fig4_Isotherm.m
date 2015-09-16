function ThesisNanoscale_Fig4_Isotherm()

    AddPaths('ThesisNanoscale');

    PhysArea = struct('N',[1,250],...
                      'L1',4,'L2',2,...
                      'alpha_deg',90);

    V2Num   = struct('Fex','SplitAnnulus','N',[80,80]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'V2Num',V2Num,...
                     'maxComp_y2',-1,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.856);%1.375);%1.25)s;
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);     

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.03,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    ComputeIsotherm(config);
    
    
    %For 45 degrees CA
    config.optsPhys.V1.epsilon_w = 1.155;
    ComputeIsotherm(config);
    
    %For 135 degrees CA
    config.optsPhys.V1.epsilon_w = 0.453;
    ComputeIsotherm(config);
    
    
    function ComputeIsotherm(config)
        CL = ContactLineHS(config);
        CL.Preprocess();    
        CL.ComputeAdsorptionIsotherm(600,'drying');    %wetting    
        CL.FittingAdsorptionIsotherm([10 14],1)
        if(config.optsPhys.kBT == 0.75)
            CL.SumRule_AdsorptionIsotherm(0.2853);
        end
    end
 
end