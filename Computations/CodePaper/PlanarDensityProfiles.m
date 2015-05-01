  function PlanarDensityProfiles()   
    
    Phys_Area = struct('shape','HalfSpace_FMT',...
                       'N',[1;60],'L1',4,'L2',2.,'y2wall',0.,...
                       'N2bound',16,'h',1,'L2_AD',2.,...
                       'alpha',pi/2);

    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',6,'N2',100);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',34,'N2disc',34);
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num);

    V1       = struct('V1DV1','zeroPotential');                  

    optsPhys = struct('V1',V1,...
                      'kBT',1,...
                      'eta',0.3744,...%0.4257,...%0.4783
                      'sigmaS',1,... 
                      'nSpecies',1);

    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;

    AddPaths();    
    close all;  
    
    config = v2struct(optsPhys,optsNum);
         
%     EX     = DDFT_2D(config);
%     EX.Preprocess();
%     CheckAverageDensities_Rosenfeld_3D(EX.IDC,EX.IntMatrFex);
% 	optsPhys.rho_iguess = optsPhys.eta*6/pi;
%     %FMT_1D_HardWall(EX.IDC,EX.IntMatrFex,optsPhys,optsNum);
%     FMT_1D(EX.IDC,EX.IntMatrFex,optsPhys,optsNum.FexNum,[],true);
    
    %*******************************
    
    optsNum.PhysArea.alpha_deg = 90;
    optsNum.V2Num = struct('Fex','SplitAnnulus','N',[80,80]);
    optsPhys.V2   = struct('V2DV2','BarkerHendersonCutoff_2D','epsilon',1,'LJsigma',1,'r_cutoff',2.5);
    
    
    optsPhys.V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',0.95);%1.375);%1.25)s;    
    optsPhys.kBT = 0.75;
    optsPhys.Dmu = 0.0;
    optsPhys = rmfield(optsPhys,'eta');
    config = v2struct(optsPhys,optsNum);
    
	CL     = ContactLineHS(config);
    CL.Preprocess();    
    CL.Compute1D('WL');
    CL.Compute1D('WG');
    %optsPhys.rho_iguess = optsPhys.eta*6/pi;
    %FMT_1D_HardWall(EX.IDC,EX.IntMatrFex,optsPhys,optsNum);
   % FMT_1D(EX.IDC,EX.IntMatrFex,EX.optsPhys,EX.optsNum.FexNum,[],true);
        
    
end