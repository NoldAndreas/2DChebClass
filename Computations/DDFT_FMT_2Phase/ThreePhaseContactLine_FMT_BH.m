function CL = ThreePhaseContactLine_FMT_BH(configIn)
    AddPaths();
    global dirData
    
    if((nargin > 0) && islogical(configIn))
        PhysArea = struct('N',[40,40],'L1_Skewed',2,'L2',2,'y2wall',0.,...
                          'N2bound',24,'h',1,'L2_AD',2.,'alpha_deg',90);

        PhysArea.Conv  = struct('L',1,'L2',[],'N',[20,20]);

        Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                           'Ncircle',1,'N1disc',35,'N2disc',34);

        optsNum = struct('PhysArea',PhysArea,...
                         'FexNum',Fex_Num,...
                         'maxComp_y2',10,...
                         'y1Shift',0,...c 
                         'plotTimes_T',100,...
                         'plotTimes_Interval',0.1);

        V1 = struct('V1DV1','Vext_Cart_7',...
                            'epsilon_w',0.74,... %0.482218,...
                            'epsilon_w_max',0.3,....
                            'tau',20);
                        
        V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'R',1);

        optsPhys = struct('V1',V1,'V2',V2,...                      
                          'kBT',0.75,...   
                          'gamma',2,...
                          'Inertial',1,...
                          'Dmu',0.0,'nSpecies',1,...
                          'sigmaS',1);      
                      
        configuration = v2struct(optsNum,optsPhys);
        CL  = ContactLine(configuration);
    elseif(nargin == 0)
        ChangeDirData([dirData filesep 'FMT_CLEq_BH_40X40_epw']);
        CL  = ContactLine();
    else
        CL  = ContactLine(configIn);
    end

    CL.GotoSubDir();    
    CL.Compute();
    
end