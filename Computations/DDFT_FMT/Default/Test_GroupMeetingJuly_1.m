function [optsNum,optsPhys,optsPlot,name] = Test_GroupMeetingJuly_1()

    Phys_Area = struct('shape','HalfSpace_FMT','N',[40,40],...
                        'L1',2,'L2',2,'L2_AD',2,...
                        'y2wall',0.,'N2bound',10,'h',1);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',6,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'plotTimes',0:0.1:2); 
                 %'Tmax',0.78,'TN',50,...
                     %'name','default'); 
                 %'DDFTCode','DDFT_DiffusionHalfSpace_NSpecies',...
                     
    %FexMatrices = {'Polar_SpectralFourier_FMTMatricesFull', ...
%               'Polar_SpectralFourier_FMTMatricesFull_Roth', ...
%               'Polar_SpectralFourier_FMTMatricesFull_3D'};
%    Fex         = {'Polar_SpectralFourier_FMT', ...
%               'Polar_SpectralFourier_FMT_Roth', ...
%               'Polar_SpectralFourier_FMT_3D'};                                  

    sigmaS = 1;%   1.1 ;
              %1.1 1.2 ];        
              
    V1       = struct('V1DV1','Vext_Cart_6','V0',1,'Axis','y2',...
                      'grav_start',-0.2,'grav_end',0,'T',5,...
                      'L1',1  ,'L2',1);                                          
    
    optsPhys = struct('V1',V1,'kBT',1,...  
                      'sigmaS',sigmaS); %[10;10] 
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
    
    eta         = 0.4257;
    fBulk       = str2func(['FexBulk_',optsNum.FexNum.Fex]);
	rhoBulk     = eta*6/pi;
    optsPhys.mu = optsPhys.kBT*log(rhoBulk) + fBulk(rhoBulk,optsPhys.kBT);            
    
    AddPaths();
    EX     = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    EX.ComputeEquilibrium();   
    EX.ComputeDynamics();
    
    %PlotRosenfeldFMT_AverageDensitiesInf(HS,IntMatrFex(1),rho_ic1D);

end                 

