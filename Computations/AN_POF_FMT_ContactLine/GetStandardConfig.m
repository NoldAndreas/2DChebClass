function config = GetStandardConfig(opts)

    global dirData
    AddPaths();    
    
    if(nargin == 0)
        alpha_deg = 60;
        epw       = 1.25;
        bounds1   = [-10 10];
        bounds2   = [0.5 18];                
        maxComp_y2 = 35;
        N          = [45,75];
    else
        v2struct(opts);
    end    

    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');        

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