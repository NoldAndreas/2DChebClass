  function [EX,res] = CheckFMTComputation()

    disp('** CheckFMTComputationSkewed **');
    
    Phys_Area = struct('shape','HalfSpace_FMT',...
                       'N',[1;70],'L1',2,'L2',2.,'y2wall',0.,...
                       'N2bound',30,'h',1,'L2_AD',2.,...
                       'alpha',pi/2);

    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',6,'N2',100);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num);

    V1       = struct('V1DV1','zeroPotential');                  

    optsPhys = struct('V1',V1,...
                      'kBT',1,...
                      'eta',0.4257,...%0.4783
                      'sigmaS',1,... 
                      'nSpecies',1);

    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;

    AddPaths();    
    close all;  
         
    EX   = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    CheckAverageDensities_Rosenfeld_3D(EX.IDC,EX.IntMatrFex);
	optsPhys.rho_iguess = optsPhys.eta*6/pi;
    FMT_1D_HardWall(EX.IDC,EX.IntMatrFex,optsPhys,optsNum);
    res.fig_handles{1} = gcf; 
end