 function config = ThesisNanoscale_GetStandardConfig(alpha_deg,epw)
 
    if(nargin < 1)
        alpha_deg = [];
    end
    if(nargin < 2)
        epw = 0;
    end
            
    PhysArea = struct('N',[45 75],'L1',4,'L2',2,...
                      'alpha_deg',alpha_deg);
    
    V2Num     = struct('Fex','SplitAnnulus','N',[50,50]);
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid','Ncircle',1,'N1disc',30,'N2disc',30); 

    optsNum = struct('PhysArea',PhysArea,...%'PlotAreaCart',PlotArea,...                     
                     'FexNum',Fex_Num,...
                     'V2Num',V2Num,...
                     'maxComp_y2',25);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',epw);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        

end