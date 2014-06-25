function PlotExemplaricGrid()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'POF_FMT_ContactLine'],'ORG');    

    PhysArea = struct('N',[45,75],...
                      'L1',4,'L2',2,'L2_AD',2.,...
                      'y2wall',0.,...
                      'N2bound',14,'h',1,...
                      'alpha_deg',40);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);
    
    PlotArea = struct('y1Min',-5,'y1Max',10,...
                      'y2Min',0.5,'y2Max',10,'N1',100,'N2',100);    

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50); %35,34
 
    optsNum = struct('PhysArea',PhysArea,'FexNum',Fex_Num,...
                     'PlotArea',PlotArea,'maxComp_y2',35,'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.375);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %***********************************************************        
    filename   = [dirData,filesep,'Job40_DisjoiningPressure_',getTimeStr(),'.txt'];       
    Struct2File(filename,config,['Computed at ',datestr(now)]);
        
  %  ConvergenceSurfaceTensions(config);
        
    close all;
    ChangeDirData();                
	
    CLT = ContactLine(config);     
    CLT.Preprocess();
    CLT.ComputeEquilibrium();  
    
    f1 = figure('color','white','Position',[0 0 800 500]);            
        
    CLT.InitInterpolation(true);
        
	optDetails.y2CartShift = -0.5;
    optDetails.clabel      = false;  
    optDetails.linewidth   = 2;  

    rhoLiq_sat    = CLT.optsPhys.rhoLiq_sat;
	rhoGas_sat    = CLT.optsPhys.rhoGas_sat;    
    drho                 = rhoLiq_sat - rhoGas_sat;
    optDetails.nContours = rhoGas_sat + 0.5*drho;
    optDetails.linecolor = 'k';    
  
    CLT.HS.doPlots(CLT.rho_eq,'contour',optDetails);  hold on;  
    
    optDetails.nthGridLines = 2;
    CLT.HS.PlotGridLines(optDetails); 
    
    xlabel('$x/\sigma$','Interpreter','Latex','fontsize',25);
	ylabel('$y/\sigma$','Interpreter','Latex','fontsize',25);                
    
	print2eps([dirData filesep 'ExemplaricGrid'],f1);
    saveas(f1,[dirData filesep 'ExemplaricGrid.fig']);

end