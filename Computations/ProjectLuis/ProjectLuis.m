%function ProjectLuis
   global dirData
   AddPaths();        
   ChangeDirData([dirData filesep 'MMNP'],'ORG');    
   close all;
                  
   bounds    = [-10,10]; %[0,20]
   %alpha_deg = 120; epw      = 0.7; %40 %1.375     
   alpha_deg = 135; epw      = 0.55; 
               
    PhysArea = struct('N',[50 80],'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',alpha_deg);

    PlotAreaCart = struct('y1Min',bounds(1),'y1Max',bounds(2),...
                      'y2Min',0.5,'y2Max',15.5,...
                      'zMax',4,...
                      'N1',100,'N2',100);

    V2Num   = struct('Fex','SplitDisk','L',1,'L2',[],'N',[34,34]);    
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50); %35,34

    optsNum = struct('PhysArea',PhysArea,...
                     'PlotAreaCart',PlotAreaCart,...
                     'FexNum',Fex_Num,...
                     'V2Num',V2Num,...
                     'maxComp_y2',35);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',epw);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);     
    
    
    shapeSL = struct('yMin',PlotAreaCart.y1Min,...
                     'yMax',PlotAreaCart.y1Max,...
                     'N',150);                           
    %*****************************************************    

    CLT = ContactLineHS(config);     
    CLT.Preprocess();        
    CLT.ComputeEquilibrium();      
    
    CLT.y1_SpectralLine = SpectralLine(shapeSL);
    CLT.y1_SpectralLine.ComputeAll();

    CLT.Compute_hContour(0.5);
    CLT.Compute_hIII();
    CLT.PlotDensitySlicesNormalInterface();
    %CLT.PlotDensitySlices();
    %CLT.PlotDisjoiningPressures();      
   

%end