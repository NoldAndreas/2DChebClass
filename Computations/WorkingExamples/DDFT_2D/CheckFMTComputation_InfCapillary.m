  function [EX,res] = CheckFMTComputation_InfCapillary()

    disp('** CheckFMTComputationSkewed **');
    
    Phys_Area = struct('shape','InfCapillary_FMT',...
                       'N',[1;50],...
                       'L1',2,...
                       'y2Min',0.5,'y2Max',8,...
                       'N2bound',30);

    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',Phys_Area.y2Max,'N2',100);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',30,'N2disc',30);
    
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
    CheckAverageDensities_Rosenfeld_3D(EX.IDC,EX.IntMatrFex,'checkTop');
	optsPhys.rho_iguess = optsPhys.eta*6/pi;
    %FMT_1D_HardWall(EX.IDC,EX.IntMatrFex,optsPhys,optsNum);
    res.fig_handles{1} = gcf;
end